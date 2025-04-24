import argparse
import os, yaml, torch, wandb
import numpy as np
import pandas as pd

from utils.logging import get_logger
from utils.dataset import get_loader
from model import get_model
from utils.inference import inference_epoch

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
logger = get_logger(__name__)

parser = argparse.ArgumentParser()
parser.add_argument('--model_dir', type=str, required=True)
parser.add_argument('--ckpt', type=str, required=True)
parser.add_argument('--splits', type=str, required=True)
parser.add_argument('--split_key', type=str, default=None)
parser.add_argument('--inf_mols', type=int, default=1000)
parser.add_argument('--wandb', type=str, default=None)
parser.add_argument('--num_workers', type=int, default=None)
parser.add_argument('--ode', action='store_true', default=False)
parser.add_argument('--elbo', action='store_true', default=False)
parser.add_argument('--alpha', type=float, default=0)
parser.add_argument('--beta', type=float, default=1)
parser.add_argument('--num_samples', type=int, default=1)
parser.add_argument('--inf_step', type=float, default=0.5)
parser.add_argument('--elbo_step', type=float, default=0.2)
parser.add_argument('--inf_type', type=str, choices=['entropy', 'rate'], default='rate')
parser.add_argument('--max_len', type=int, default=1024)
parser.add_argument('--inf_Hf', type=float, default=None)
parser.add_argument('--inf_kmin', type=int, default=None)
parser.add_argument('--inf_tmin', type=float, default=None)
parser.add_argument('--inf_cutoff', type=int, default=None)
parser.add_argument('--embeddings_dir', type=str, default=None)
parser.add_argument('--pdb_dir', type=str, default=None)
parser.add_argument('--embeddings_key', type=str, default=None, choices=['name', 'reference'])

inf_args = parser.parse_args()

with open(f'{inf_args.model_dir}/args.yaml') as f:
    args = argparse.Namespace(**yaml.full_load(f))

args.wandb = inf_args.wandb
if inf_args.wandb:
    wandb.login(key=os.environ['WANDB_API_KEY'])
    wandb.init(
        entity=os.environ['WANDB_ENTITY'],
        settings=wandb.Settings(start_method="fork"),
        project=args.wandb,
        name=str(args.time),
        config=args.__dict__ | inf_args.__dict__
    )

# Overwrite with inference-time args
args.splits = inf_args.splits
args.inf_mols = inf_args.inf_mols
args.inf_type = inf_args.inf_type
args.num_samples = inf_args.num_samples
args.inf_step = inf_args.inf_step
args.max_len = inf_args.max_len
args.ode = inf_args.ode
args.alpha, args.beta = inf_args.alpha, inf_args.beta
args.elbo_step = inf_args.elbo_step
args.inference_mode = True

if inf_args.num_workers is not None:
    args.num_workers = inf_args.num_workers
if inf_args.inf_Hf: args.train_Hf = inf_args.inf_Hf
if inf_args.inf_kmin: args.train_kmin = inf_args.inf_kmin
if inf_args.inf_tmin: args.train_tmin = inf_args.inf_tmin
if inf_args.pdb_dir: args.pdb_dir = inf_args.pdb_dir
if inf_args.embeddings_dir: args.embeddings_dir = inf_args.embeddings_dir
if inf_args.embeddings_key: args.embeddings_key = inf_args.embeddings_key

def main():
    logger.info(f'Loading splits {args.splits}')
    try:
        splits = pd.read_csv(args.splits).set_index('path')
    except:
        splits = pd.read_csv(args.splits).set_index('name')

    logger.info("Constructing model")
    model = get_model(args).to(device)
    ckpt = os.path.join(inf_args.model_dir, inf_args.ckpt)

    logger.info(f'Loading weights from {ckpt}')
    state_dict = torch.load(ckpt, map_location=torch.device('cpu'))
    model.load_state_dict(state_dict['model'], strict=True)
    ep = state_dict['epoch']

    val_loader = get_loader(args, None, splits, mode=inf_args.split_key, shuffle=False)
    samples, log = inference_epoch(args, model, val_loader.dataset, device=device, pdbs=True, elbo=inf_args.elbo)

    means = {key: np.mean(log[key]) for key in log if key != 'path'}
    logger.info(f"Inference epoch {ep}: len {len(log.get('rmsd', []))} MEANS {means}")

    # Safe and unique filename
    split_name = os.path.basename(args.splits).replace('.csv', '')
    inf_name = f"{split_name}.ep{ep}.num{args.num_samples}.step{args.inf_step}.alpha{args.alpha}.beta{args.beta}"
    if inf_args.inf_step != inf_args.elbo_step:
        inf_name += f".elbo{args.elbo_step}"

    # Save results
    output_dir = os.path.join(os.getcwd(), "outputs")
    os.makedirs(output_dir, exist_ok=True)
    csv_path = os.path.join(output_dir, f"{inf_name}.csv")
    pd.DataFrame(log).set_index('path').to_csv(csv_path)
    logger.info(f"Saved inference results to {csv_path}")

    # Save PDBs
    save_dir = os.path.join(inf_args.model_dir, inf_name)
    os.makedirs(save_dir, exist_ok=True)
    for samp in samples:
        samp.pdb.write(os.path.join(save_dir, f"{os.path.basename(samp.path)}.{samp.copy}.anim.pdb"), reverse=True)
        samp.pdb.clear().add(samp.Y).write(os.path.join(save_dir, f"{os.path.basename(samp.path)}.{samp.copy}.pdb"))

if __name__ == '__main__':
    main()

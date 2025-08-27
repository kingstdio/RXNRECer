import sys,os
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, f'{project_root}/')
from rxnrecer.config import config as cfg
from types import SimpleNamespace
from rxnrecer.lib.model import mactive as Mactive





def train(args):
    
    Mactive.train(
        train_path=args.train_path, 
        test_path=args.test_path,
        batch_size=args.batch_size, 
        esm_out_dim=args.esm_out_dim, 
        gru_h_dim=args.gru_h_dim, 
        att_dim=args.att_dim, 
        dropout_rate=args.dropout_rate, 
        freeze_esm_layers = args.freeze_esm_layers,
        epoches=args.epoches,
        args=args
    )
    


if __name__ == '__main__':
    

    args_dict = {
        'train_path': cfg.FILE_DS_TRAIN,
        'test_path': cfg.FILE_DS_TEST,
        'batch_size': 8,
        'esm_out_dim': 1280,
        'gru_h_dim': 256,
        'att_dim': 32,
        'dropout_rate': 0.2,
        'freeze_esm_layers': 23,
        'epoches': 5,
        'taskname': 'train_sample',
        'devices': 2,
    }

    args = SimpleNamespace(**args_dict)

    train(args)
    print('ok')
import sys
import os
import argparse
import configparser
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from POSTOA_main import POSTOA_main


def run(year, month, input_dir, output_dir, config_fname=None, varname=None):
    cfg = configparser.ConfigParser(interpolation=None)
    cfg.read(Path(__file__).parent / 'config.ini')

    if varname is None:
        varname = cfg['analysis']['var_name']
    if config_fname is None:
        config_fname = cfg['paths']['config_fname']

    month_str = f'{month:02d}'
    climpath = cfg['paths']['climpath']

    clim = cfg['paths']['clim_template'].format(
        climpath=climpath, year=str(year), month=month_str, var=varname)
    nam_atlas_pres = cfg['paths']['atlas_pres_template'].format(
        climpath=climpath, month=month_str)

    POSTOA_main(config_fname, [15, month, year], varname, 0,
                nam_atlas_pres, clim,
                input_dir=input_dir, output_dir=output_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='POSTOA launcher — post-process OA results for ISAS')
    parser.add_argument('--input',  required=True,
                        help='Directory containing PREOA zone files')
    parser.add_argument('--output', required=True,
                        help='Root directory for POSTOA output')
    parser.add_argument('--year',   required=True, type=int,
                        help='Analysis year (e.g. 2010)')
    parser.add_argument('--month',  required=True, type=int,
                        help='Analysis month (1-12, or zero-padded e.g. 01)')
    parser.add_argument('--config', default=None,
                        help='XML configuration file (default: from config.ini)')
    parser.add_argument('--var',    default=None,
                        help='Variable name to process (default: from config.ini)')
    args = parser.parse_args()

    run(
        year=args.year,
        month=args.month,
        input_dir=args.input,
        output_dir=args.output,
        config_fname=args.config,
        varname=args.var,
    )

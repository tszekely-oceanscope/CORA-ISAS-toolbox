#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import configparser
import logging
import time
import os
import sys
from pathlib import Path
from collections import defaultdict


def write_config(isasconfig):
    """
    Write the .cnf config file from the .ini config.

    Parameters
    ----------
    isasconfig : ConfigParser

    Returns
    -------
    cfgfilename : str - path to the generated .cnf file
    """
    tmp_dir = isasconfig['path']['tmp']
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)

    cfgfilename = f"{tmp_dir}/tmp-{time.strftime('%Y%m%d-%H%M%S')}.cnf"
    print(f'ConfigGen : writing temp configfile @ {cfgfilename}')

    with open(cfgfilename, 'w') as f:
        f.write(f"{isasconfig['config']['param']}\n")
        f.write(f"{isasconfig['config']['clim']}\n")
        f.write(f"{isasconfig['config']['cov']}\n")
        f.write(f"{isasconfig['config']['bathy']}\n")
        f.write("{} {} {}\n".format(
            isasconfig['config']['covar_ls_x'],
            isasconfig['config']['covar_ls_y'],
            isasconfig['config']['covar_ls_t']))
        f.write(f"{isasconfig['config']['covar_ms_t']}\n")
        f.write("{} {} {}\n".format(
            isasconfig['config']['var_weigh_ls'],
            isasconfig['config']['var_weigh_ms'],
            isasconfig['config']['var_weigh_ur']))
        f.write("{} {} {} {}\n".format(
            isasconfig['config']['cd_x'],
            isasconfig['config']['cd_y'],
            isasconfig['config']['cd_z'],
            isasconfig['config']['cd_t']))
        f.write(f"{isasconfig['config']['variance']}\n")
        f.write("{} {}\n".format(
            isasconfig['config']['qc_max'],
            isasconfig['config']['mx_std']))
        f.write(f"{isasconfig['config']['cov_max']}\n")

    return cfgfilename


def define_tasks(isasconfig):
    """
    Build the list of valid tasks from the series config.

    Iterates over years, months, bassins, and boxes defined in [series],
    skips boxes not in real_box, and skips input files that don't exist on disk.

    Parameters
    ----------
    isasconfig : ConfigParser

    Returns
    -------
    tasks : list of [y, m, box, inputfile, job_data_dir]
    """
    series = isasconfig['series']
    param  = isasconfig['config']['param']

    real_box = (list(range(101, 162)) +
                list(range(201, 239)) +
                list(range(301, 383)) +
                list(range(401, 405)) +
                list(range(501, 507)))

    base_data_dir = isasconfig['path']['data']

    tasks = []
    for y in range(int(series['first_y']), int(series['last_y']) + 1):
        for m in range(int(series['first_m']), int(series['last_m']) + 1):
            mstr = f'{m:02d}'
            job_data_dir = os.path.join(base_data_dir, param, str(y), mstr)
            for b in range(int(series['first_bassin']), int(series['last_bassin']) + 1):
                first_box = int(f"{b}{series['first_box']}")
                last_box  = int(f"{b}{series['last_box']}")
                for box in range(first_box, last_box + 1):
                    if box not in real_box:
                        continue
                    inputfile = f"OA_{y}{mstr}15_{box}_dat_{param}.nc"
                    fpath = os.path.join(job_data_dir, inputfile)
                    if not os.path.isfile(fpath):
                        logging.warning(f'DefineTasks: SKIP {inputfile} — not found in {job_data_dir}')
                    else:
                        if [y, m, box, inputfile, job_data_dir] not in tasks:
                            tasks.append([y, m, box, inputfile, job_data_dir])

    print(f'Define_tasks: {len(tasks)} valid tasks total')
    return tasks


def write_infile(tasks, isasconfig, cfgfilename, ts):
    """
    Write one .in file per data directory (grouped by job_data_dir).

    Each .in file contains the data dir, config path, task count,
    and the list of input NetCDF filenames for that directory.

    Parameters
    ----------
    tasks      : list of [y, m, box, inputfile, job_data_dir]
    isasconfig : ConfigParser
    cfgfilename: str - path to the .cnf config file
    ts         : str - timestamp string

    Returns
    -------
    infiles : list of str - paths to the generated .in files
    """
    tmp_dir = isasconfig['path']['tmp']
    Path('err').mkdir(exist_ok=True)

    groups = defaultdict(list)
    for t in tasks:
        groups[t[4]].append(t)

    infiles = []
    for i, (data_dir, group_tasks) in enumerate(groups.items()):
        infilename = os.path.join(tmp_dir, f'job-{ts}-{i:03d}.in')
        with open(infilename, 'w') as f:
            f.write(f"{data_dir}/$\n")
            f.write(f"{data_dir}/$\n")
            f.write(f"{cfgfilename}$\n")
            f.write(f"{len(group_tasks)}\n")
            for t in group_tasks:
                f.write(f"{t[3]}\n")
        infiles.append(infilename)
        print(f'  .in: {infilename} ({len(group_tasks)} tasks, dir: {data_dir})')

    return infiles


def main(mode):
    """
    Entry point. Reads config, builds tasks, writes .in files, and produces list.lst.

    Parameters
    ----------
    mode : str - 'TEMP' or 'PSAL'
    """
    config_map = {
        'TEMP': ('model_config_sciT.ini', './log/runpy_co_sciT.log'),
        'PSAL': ('model_config_sciS.ini', './log/runpy_co_sciS.log'),
    }
    if mode not in config_map:
        print(f'Unknown mode: {mode}. Use TEMP or PSAL.')
        sys.exit(1)

    configname, logname = config_map[mode]
    Path('./log').mkdir(exist_ok=True)
    logging.basicConfig(filename=logname, level=logging.INFO)

    isasconfig = configparser.ConfigParser()
    try:
        isasconfig.read(configname)
    except Exception as e:
        logging.warning(f'Error reading config file: {e}')
        sys.exit(1)

    print("data dir:", isasconfig['path']['data'])
    print("param:", isasconfig['config']['param'])

    cfgfilename = write_config(isasconfig)
    tasks       = define_tasks(isasconfig)

    if not tasks:
        print('No tasks found.')
        return

    ts      = time.strftime('%Y%m%d-%H%M%S')
    infiles = write_infile(tasks, isasconfig, cfgfilename, ts)

    with open('./list.lst', 'w') as f:
        for infilename in infiles:
            f.write(f"../isas_f90/OA_main ./{infilename}\n")

    print(f'list.lst written ({len(infiles)} entries)')


main('PSAL')

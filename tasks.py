#!/usr/bin/env python

import json
import os
import re
import sys

from doit.tools import run_once, create_folder, title_with_actions
from doit.task import clean_targets, dict_to_task

import screed
import khmer

def clean_folder(target):
    try:
        rmtree(target)
    except OSError:
        pass

seq_ext = re.compile(r'(.fasta)|(.fa)|(.fastq)|(.fq)')
def strip_seq_extension(fn):
    return seq_ext.split(fn)[0]

def create_task_object(task_dict_func):
    '''Wrapper to decorate functions returning pydoit
    Task dictionaries and have them return pydoit Task
    objects
    '''
    def d_to_t(*args, **kwargs):
        ret_dict = task_dict_func(*args, **kwargs)
        return dict_to_task(ret_dict)
    return d_to_t

@create_task_object
def get_sample_randomly_task(sample_fn, target_fn, n_reads):

    cmd = 'sample-reads-randomly.py -N {n_reads} '\
          '-o {target_fn} {sample_fn}'.format(**locals())

    return {'name': 'random-sample:' + os.path.basename(sample_fn),
            'actions': [cmd],
            'file_dep': [sample_fn],
            'targets': [target_fn],
            'clean': [clean_targets]}

@create_task_object
def get_gzip_task(fn):

    cmd = 'gzip {fn}'.format(fn=fn)

    return {'name': 'gzip:' + os.path.basename(fn),
            'actions': [cmd],
            'file_dep': [fn],
            'targets': [fn + '.gz'],
            'clean': [clean_targets]}

@create_task_object
def download_task(url, target_fn, label='default'):

    cmd = 'curl -o {target_fn} {url}'.format(**locals())
    name = 'download_gunzip:' + target_fn

    return {'title': title_with_actions,
            'name': name,
            'actions': [cmd],
            'targets': [target_fn],
            'clean': [clean_targets],
            'uptodate': [run_once]}

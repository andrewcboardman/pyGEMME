from typing import List, Tuple, Dict, Any, Optional, Union
import requests
import logging
import time
import os
from tqdm import tqdm
import random
import tarfile
from pathlib import Path
from io import StringIO
from Bio import SeqIO, AlignIO
import re

logger = logging.getLogger(__name__)

DEFAULT_API_SERVER = "https://api.colabfold.com"

ACCEPT_DEFAULT_TERMS = \
"""
WARNING: You are welcome to use the default MSA server, however keep in mind that it's a
limited shared resource only capable of processing a few thousand MSAs per day. Please
submit jobs only from a single IP address. We reserve the right to limit access to the
server case-by-case when usage exceeds fair use. If you require more MSAs: You can 
precompute all MSAs with `colabfold_search` or host your own API and pass it to `--host-url`
"""

TQDM_BAR_FORMAT = '{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'


def get_query_sequence_uniprot(uniprot_name, results_dir):
  if isinstance(results_dir, str):
    results_dir = Path(results_dir)
  
  url = f'https://rest.uniprot.org/uniprotkb/{uniprot_name}.fasta'
  r = requests.get(url)
  if r.ok:
    fasta_io = StringIO(r.text) 
    seq_record = SeqIO.read(fasta_io, "fasta") 
    SeqIO.write(seq_record, results_dir.joinpath(f'{uniprot_name}.fasta'),format='fasta')
    return str(seq_record.seq)
  else:
    raise Exception('Error retrieving target sequence from Uniprot!')
  
def remove_inserts_a3m(query_name, raw_msa_file):
  
  raw_msa = open(raw_msa_file,'r').readlines()
  clean_msa = []
  for i, line in enumerate(raw_msa):
    if i == 0:
      clean_msa.append('>' + query_name)
    elif line[0] == '>':
      header = line.split('\t')
      name = f'{header[0]}/{header[-3]}-{header[-2]}'
      clean_msa.append(name)
    else:
      aligned_seq_without_inserts = re.sub(r'[^A-Z_.-]',r'',line)
      clean_msa.append(aligned_seq_without_inserts)
  return [line + '\n' for line in clean_msa]


def run_mmseqs2(x, prefix, use_env=True, use_filter=True,
                use_templates=False, filter=None, use_pairing=False, pairing_strategy="greedy",
                host_url="https://api.colabfold.com") -> Tuple[List[str], List[str]]:
  
  submission_endpoint = "ticket/pair" if use_pairing else "ticket/msa"

  def submit(seqs, mode, N=101):
    n, query = N, ""
    for seq in seqs:
      query += f">{n}\n{seq}\n"
      n += 1

    while True:
      error_count = 0
      try:
        # https://requests.readthedocs.io/en/latest/user/advanced/#advanced
        # "good practice to set connect timeouts to slightly larger than a multiple of 3"
        res = requests.post(f'{host_url}/{submission_endpoint}', data={'q':query,'mode': mode}, timeout=6.02)
      except requests.exceptions.Timeout:
        logger.warning("Timeout while submitting to MSA server. Retrying...")
        continue
      except Exception as e:
        error_count += 1
        logger.warning(f"Error while fetching result from MSA server. Retrying... ({error_count}/5)")
        logger.warning(f"Error: {e}")
        time.sleep(5)
        if error_count > 5:
          raise
        continue
      break

    try:
      out = res.json()
    except ValueError:
      logger.error(f"Server didn't reply with json: {res.text}")
      out = {"status":"ERROR"}
    return out

  def status(ID):
    while True:
      error_count = 0
      try:
        res = requests.get(f'{host_url}/ticket/{ID}', timeout=6.02)
      except requests.exceptions.Timeout:
        logger.warning("Timeout while fetching status from MSA server. Retrying...")
        continue
      except Exception as e:
        error_count += 1
        logger.warning(f"Error while fetching result from MSA server. Retrying... ({error_count}/5)")
        logger.warning(f"Error: {e}")
        time.sleep(5)
        if error_count > 5:
          raise
        continue
      break
    try:
      out = res.json()
    except ValueError:
      logger.error(f"Server didn't reply with json: {res.text}")
      out = {"status":"ERROR"}
    return out

  def download(ID, path):
    error_count = 0
    while True:
      try:
        res = requests.get(f'{host_url}/result/download/{ID}', timeout=6.02)
      except requests.exceptions.Timeout:
        logger.warning("Timeout while fetching result from MSA server. Retrying...")
        continue
      except Exception as e:
        error_count += 1
        logger.warning(f"Error while fetching result from MSA server. Retrying... ({error_count}/5)")
        logger.warning(f"Error: {e}")
        time.sleep(5)
        if error_count > 5:
          raise
        continue
      break
    with open(path,"wb") as out: out.write(res.content)

  # process input x
  seqs = [x] if isinstance(x, str) else x

  # compatibility to old option
  if filter is not None:
    use_filter = filter

  # setup mode
  if use_filter:
    mode = "env" if use_env else "all"
  else:
    mode = "env-nofilter" if use_env else "nofilter"

  if use_pairing:
    use_templates = False
    use_env = False
    mode = ""
    # greedy is default, complete was the previous behavior
    if pairing_strategy == "greedy":
      mode = "pairgreedy"
    elif pairing_strategy == "complete":
      mode = "paircomplete"

  # define path
  path = f"{prefix}_{mode}"
  if not os.path.isdir(path): os.mkdir(path)

  # call mmseqs2 api
  tar_gz_file = f'{path}/out.tar.gz'
  N,REDO = 101,True

  # deduplicate and keep track of order
  seqs_unique = []
  #TODO this might be slow for large sets
  [seqs_unique.append(x) for x in seqs if x not in seqs_unique]
  Ms = [N + seqs_unique.index(seq) for seq in seqs]
  # lets do it!
  if not os.path.isfile(tar_gz_file):
    TIME_ESTIMATE = 150 * len(seqs_unique)
    with tqdm(total=TIME_ESTIMATE, bar_format=TQDM_BAR_FORMAT) as pbar:
      while REDO:
        pbar.set_description("SUBMIT")

        # Resubmit job until it goes through
        out = submit(seqs_unique, mode, N)
        while out["status"] in ["UNKNOWN", "RATELIMIT"]:
          sleep_time = 5 + random.randint(0, 5)
          logger.error(f"Sleeping for {sleep_time}s. Reason: {out['status']}")
          # resubmit
          time.sleep(sleep_time)
          out = submit(seqs_unique, mode, N)

        if out["status"] == "ERROR":
          raise Exception(f'MMseqs2 API is giving errors. Please confirm your input is a valid protein sequence. If error persists, please try again an hour later.')

        if out["status"] == "MAINTENANCE":
          raise Exception(f'MMseqs2 API is undergoing maintenance. Please try again in a few minutes.')

        # wait for job to finish
        ID,TIME = out["id"],0
        pbar.set_description(out["status"])
        while out["status"] in ["UNKNOWN","RUNNING","PENDING"]:
          t = 5 + random.randint(0,5)
          logger.error(f"Sleeping for {t}s. Reason: {out['status']}")
          time.sleep(t)
          out = status(ID)
          pbar.set_description(out["status"])
          if out["status"] == "RUNNING":
            TIME += t
            pbar.update(n=t)
          #if TIME > 900 and out["status"] != "COMPLETE":
          #  # something failed on the server side, need to resubmit
          #  N += 1
          #  break

        if out["status"] == "COMPLETE":
          if TIME < TIME_ESTIMATE:
            pbar.update(n=(TIME_ESTIMATE-TIME))
          REDO = False

        if out["status"] == "ERROR":
          REDO = False
          raise Exception(f'MMseqs2 API is giving errors. Please confirm your input is a valid protein sequence. If error persists, please try again an hour later.')

      # Download results
      download(ID, tar_gz_file)

  # prep list of a3m files
  if use_pairing:
    a3m_files = [f"{path}/pair.a3m"]
  else:
    a3m_files = [f"{path}/uniref.a3m"]
    if use_env: a3m_files.append(f"{path}/bfd.mgnify30.metaeuk30.smag30.a3m")

  # extract a3m files
  if any(not os.path.isfile(a3m_file) for a3m_file in a3m_files):
    with tarfile.open(tar_gz_file) as tar_gz:
      tar_gz.extractall(path)

  # templates
  if use_templates:
    templates = {}
    #print("seq\tpdb\tcid\tevalue")
    for line in open(f"{path}/pdb70.m8","r"):
      p = line.rstrip().split()
      M,pdb,qid,e_value = p[0],p[1],p[2],p[10]
      M = int(M)
      if M not in templates: templates[M] = []
      templates[M].append(pdb)
      #if len(templates[M]) <= 20:
      #  print(f"{int(M)-N}\t{pdb}\t{qid}\t{e_value}")

    template_paths = {}
    for k,TMPL in templates.items():
      TMPL_PATH = f"{prefix}_{mode}/templates_{k}"
      if not os.path.isdir(TMPL_PATH):
        os.mkdir(TMPL_PATH)
        TMPL_LINE = ",".join(TMPL[:20])
        response = None
        while True:
          error_count = 0
          try:
            # https://requests.readthedocs.io/en/latest/user/advanced/#advanced
            # "good practice to set connect timeouts to slightly larger than a multiple of 3"
            response = requests.get(f"{host_url}/template/{TMPL_LINE}", stream=True, timeout=6.02)
          except requests.exceptions.Timeout:
            logger.warning("Timeout while submitting to template server. Retrying...")
            continue
          except Exception as e:
            error_count += 1
            logger.warning(f"Error while fetching result from template server. Retrying... ({error_count}/5)")
            logger.warning(f"Error: {e}")
            time.sleep(5)
            if error_count > 5:
              raise
            continue
          break
        with tarfile.open(fileobj=response.raw, mode="r|gz") as tar:
          tar.extractall(path=TMPL_PATH)
        os.symlink("pdb70_a3m.ffindex", f"{TMPL_PATH}/pdb70_cs219.ffindex")
        with open(f"{TMPL_PATH}/pdb70_cs219.ffdata", "w") as f:
          f.write("")
      template_paths[k] = TMPL_PATH

  # gather a3m lines
  a3m_lines = {}
  for a3m_file in a3m_files:
    update_M,M = True,None
    for line in open(a3m_file,"r"):
      if len(line) > 0:
        if "\x00" in line:
          line = line.replace("\x00","")
          update_M = True
        if line.startswith(">") and update_M:
          M = int(line[1:].rstrip())
          update_M = False
          if M not in a3m_lines: a3m_lines[M] = []
        a3m_lines[M].append(line)

  # return results

  a3m_lines = ["".join(a3m_lines[n]) for n in Ms]

  if use_templates:
    template_paths_ = []
    for n in Ms:
      if n not in template_paths:
        template_paths_.append(None)
        #print(f"{n-N}\tno_templates_found")
      else:
        template_paths_.append(template_paths[n])
    template_paths = template_paths_


  return (a3m_lines, template_paths) if use_templates else a3m_lines


# file takes a fasta alignment and outputs pir and modpy files for modeller
import re
import sys
from collections import defaultdict

#
# usage
# create_modeller_inputs.py query example/example_alignment.fasta /scratch1/NOT_BACKED_UP/dbuchan/pdb/ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/ /webdata/data/cath_pdb/ 0
#


def parseQuerySeq(fasta):
    name_re = re.compile(r'^>(.+)')
    align_deets = defaultdict(str)
    seq_name = ''
    with open(fasta, 'r') as msa:
        for line in msa:
            line = line.rstrip()
            m = re.search(name_re, line)
            if m:
                m.group(1)
                if m.group(1) not in 'Query':
                    seq_name = m.group(1)
                    align_deets['KNOWN'] = m.group(1)
                else:
                    seq_name = 'Query'
            else:
                align_deets[seq_name] += line
    align_deets['Query'] = align_deets['Query'].replace("-", "")
    align_deets['Query'] = align_deets['Query'].replace("*", "")
    return(align_deets)


def printFasta(deets, job_id):
    file_name = job_id+'.fasta'
    with open(file_name, 'w') as fasta:
        fasta.write('>Query\n')
        fasta.write(deets['Query'])


def printPIR(deets, job_id):
    file_name = job_id+'.pir'
    with open(file_name, 'w') as pir:
        pir.write('>P1;Query\n')
        pir.write('sequence:Query::::::::\n')
        pir.write(deets['Query']+'*\n')
        pir.write('>P1;'+deets['KNOWN']+'\n')
        pir.write('structure:' + deets['KNOWN'][:4] + '::' + deets['KNOWN'][4:5] + '::' + deets['KNOWN'][4:5] + '::::\n')
        pir.write(deets[deets['KNOWN']]+'*\n')


def printModPY(deets, type, pdb_path, cath_dom_path, job_id):
    # type controls if we're writing for cath domains
    # 0 == pdb, 1 == cath
    file_name = job_id+'.mod.py'
    querypir = job_id+'.pir'
    with open(file_name, 'w') as modpy:
        strModPy = "from modeller import *\n"
        strModPy += "from modeller.automodel import *\n\n"
        strModPy += "env = environ()\n"

        if type == 0:
            strModPy += "env.io.atom_files_directory = \'" + pdb_path + "\'\n"
        else:
            strModPy += "env.io.atom_files_directory = \'" + cath_dom_path + "\'\n"

        strModPy += "a = automodel(env, alnfile=\'" + querypir + "\',\n"
        strModPy += "              knowns=("

        if type == 0:
            strModPy += "\'" + deets['KNOWN'][:4].lower()+deets['KNOWN'][4:6] + "\'"
        else:
            strModPy += "\'" + deets['KNOWN'][:4] + "\'"

        strModPy += "),\n"
        strModPy += "              sequence=\'Query\')\n"
        strModPy += "a.starting_model = 1\n"
        strModPy += "a.ending_model = 1\n"
        strModPy += "a.make()\n"
        strModPy += "ok_models = filter(lambda x: x[\'failure\'] is None, a.outputs)\n"
        strModPy += "key = \'molpdf\'\n"
        strModPy += "ok_models.sort(lambda a,b: cmp(a[key], b[key]))\n"
        strModPy += "m = ok_models[0]\n"
        strModPy += "print \"Top model:%s\" % m[\'name\']\n"

        modpy.write(strModPy)


job_id = sys.argv[1]  # identifier/path root for files
alignment = sys.argv[2]  # fasta alignment
pdb = sys.argv[3]  # path to pdb coordinate files
cath_pdb = sys.argv[4]  # path to cath coordinate files
model_type = sys.argv[5]  # 0 - pdb mode, 1 - cath model

queryDetails = parseQuerySeq(alignment)
printFasta(queryDetails, job_id)
printPIR(queryDetails, job_id)
printModPY(queryDetails, 0, pdb, cath_pdb, job_id)

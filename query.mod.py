from modeller import *
from modeller.automodel import *

env = environ()
env.io.atom_files_directory = '/scratch1/NOT_BACKED_UP/dbuchan/pdb/ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/'
a = automodel(env, alnfile='query.pir',
              knowns=('1m0kA0'),
              sequence='Query')
a.starting_model = 1
a.ending_model = 1
a.make()
ok_models = filter(lambda x: x['failure'] is None, a.outputs)
key = 'molpdf'
ok_models.sort(lambda a,b: cmp(a[key], b[key]))
m = ok_models[0]
print "Top model:%s" % m['name']

from __future__ import print_function

import os.path

# Build SANTA config file from the FASTA sequences and a template file.
# for information on passing parameters to custom builders, see
# http://osdir.com/ml/programming.tools.scons.user/2003-02/msg00036.html
#
def generate_SANTAConfig(target, source, env, for_signature):
    print("inside generate_BuildConfig " + str(env['GENERATION']))
    prefix = os.path.splitext(os.path.basename(str(target[0])))[0]

    # randomly select COUNT sequence from generation GENERATION of source
    cmd = 'sed \'/^\(>.*\)$/N;s/\\n/:/g\' %s | grep \'_%s_\' | shuf -n %d | sed \'s/:/\\n/g\' >sample.fa\n' % (source[0], env['GENERATION'], env['COUNT'])

    # create a new santa config file using the sampled sequences
    cmd += './mksanta.py -p %s patient_santa_template.xml sample.fa > %s' % (prefix, target[0])
    # print(cmd)
    return cmd

def generate_BEASTConfig(target, source, env, for_signature):
    prefix = os.path.splitext(os.path.basename(str(target[0])))[0]
    sources = " ".join([str(s) for s in source])
    cmd = './mkbeast.py -p %s patient_beast_template.xml %s >%s' % (prefix, sources, target[0])
    # print(cmd)
    return cmd


env = Environment(SANTAJAR = os.path.expanduser('~/src/matsen/tools/santa-sim/dist/santa.jar'))

env.Append(BUILDERS={'SantaSim': Builder(action='java -jar $SANTAJAR $SOURCE', suffix='.fa')})
env.Append(BUILDERS={'BuildSANTA': Builder(generator = generate_SANTAConfig)}, suffix='.xml')
env.Append(BUILDERS={'BuildBEAST': Builder(generator = generate_BEASTConfig)}, suffix='.xml')

SConscript('SConscript', 'env')

from __future__ import print_function

import os.path

# Build SANTA config file from the FASTA sequences and a template file.
# for information on passing parameters to custom builders, see
# http://osdir.com/ml/programming.tools.scons.user/2003-02/msg00036.html
#
def generate_SANTAConfig(target, source, env, for_signature):
    prefix = os.path.splitext(os.path.basename(str(target[0])))[0]

    # # randomly select COUNT sequence from generation GENERATION of source
    # cmd = 'sed \'/^\(>.*\)$/N;s/\\n/:/g\' %s | grep \'_%s_\' | shuf -n %d | sed \'s/:/\\n/g\' >sample.fa\n' % (source[0], env['GENERATION'], env['COUNT'])

    # create a new santa config file using the sampled sequences
    cmd = './mksanta.py -p {} patient_santa_template.xml {} {} {} > {}'.format(prefix, source[0], env['GENERATION'], env['COUNT'], target[0])
    return cmd

def generate_BEASTConfig(target, source, env, for_signature):
    sources = " ".join([str(s) for s in source])
    cmd = './mkbeast.py patient_beast_template.xml %s >%s' % (sources, target[0])
    return cmd

def generate_BEASTCommand(target, source, env, for_signature):
    # extract the basename of the target to use as the name of 
    # all the log files produced by beast.
    prefix = os.path.splitext(os.path.basename(str(target[0])))[0]
    prefix = prefix.replace("_beast","")
    cmd = 'beast -overwrite -prefix {} -beagle {}'.format(prefix, source[0])
    return cmd

def beast_targets(target, source, env):
    prefix = os.path.splitext(os.path.basename(str(source[0])))[0]
    prefix = prefix.replace("_beast","")
    prefix = prefix
    target = [ prefix+"_beast.log", prefix+".trees", prefix+".states.log"]
    return target, source


env = Environment(SANTAJAR = os.path.expanduser('~/src/matsen/tools/santa-sim/dist/santa.jar'))


env.Append(BUILDERS={'SantaSim': Builder(action='java -jar $SANTAJAR $SOURCE', suffix='.fa', src_suffix='.xml')})
env.Append(BUILDERS={'BuildSANTA': Builder(generator = generate_SANTAConfig, suffix='.xml', src_suffix='.fa')})
env.Append(BUILDERS={'BuildBEAST': Builder(generator = generate_BEASTConfig, suffix='.xml')})
env.Append(BUILDERS={'BestTree': Builder(action='treeannotator $SOURCE >$TARGET', suffix=".mcc", src_suffix='.trees')})
env.Append(BUILDERS={'Beast': Builder(generator = generate_BEASTCommand, 
			      	      suffix=".trees", src_suffix='.xml',
				      emitter = beast_targets)})

SConscript('SConscript', 'env')



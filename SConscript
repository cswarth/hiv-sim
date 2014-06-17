# Note that BuildConfig() takes parameters that describe 
# which generation to sample for seeding the next patient.
# For information on passing parameters to custom builders, see
# http://osdir.com/ml/programming.tools.scons.user/2003-02/msg00036.html
import os

Import('env')
env['ENV']['PATH'] = os.environ['PATH']


env.SantaSim('patient1', 'patient1_santa.xml')

env.BuildSANTA('patient2', 'patient1', GENERATION=200, COUNT=1)
env.SantaSim('patient2')

env.BuildSANTA('patient3', 'patient2', GENERATION=400, COUNT=1)
env.SantaSim('patient3')

env.BuildSANTA('patient4', 'patient3', GENERATION=300, COUNT=1)
env.SantaSim('patient4')

env.BuildSANTA('patient5', 'patient1', GENERATION=300, COUNT=1)
env.SantaSim('patient5')

env.BuildBEAST('combined_145_beast', ['patient1.fa', 'patient4.fa','patient5.fa'])
env.Beast('combined_145_beast')
env.BestTree('combined_145', 'combined_145_beast')

env.BuildBEAST('combined_123_beast', ['patient1.fa', 'patient2.fa','patient3.fa'])
env.Beast('combined_123_beast')
env.BestTree('combined_123', 'combined_123_beast')


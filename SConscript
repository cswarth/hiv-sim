# Note that BuildConfig() takes parameters that describe 
# which generation to sample for seeding the next patient.
# For information on passing parameters to custom builders, see
# http://osdir.com/ml/programming.tools.scons.user/2003-02/msg00036.html

Import('env')

env.SantaSim('patient1', 'patient1_santa.xml')

# env.Transmission('patient2', 'patient1.fa', GENERATION=200, COUNT=1)

env.BuildSANTA('patient2.xml', 'patient1.fa', GENERATION=200, COUNT=1)
env.SantaSim('patient2', 'patient2.xml')

env.BuildSANTA('patient3.xml', 'patient2.fa', GENERATION=400, COUNT=1)
env.SantaSim('patient3', 'patient3.xml')

env.BuildSANTA('patient4.xml', 'patient3.fa', GENERATION=300, COUNT=1)
env.SantaSim('patient4', 'patient4.xml')

env.BuildSANTA('patient5.xml', 'patient1.fa', GENERATION=300, COUNT=1)
env.SantaSim('patient5', 'patient5.xml')

env.BuildBEAST('combined_beast.xml', ['patient1.fa', 'patient4.fa','patient5.fa'])

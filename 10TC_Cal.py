#A demo for calculating the thermal conductivity of polyethylene (PE)
from radonpy.core import utils, poly 
from radonpy.ff.gaff2_mod import GAFF2_mod 
from radonpy.sim import qm 
from radonpy.sim.preset import eq, tc
import matplotlib
 
matplotlib.use('AGG')
 
smiles = '[*]CC[*]' 
ter_smiles = '*C' 
temp = 300 
press = 1.0 
omp_psi4 = 40 
mpi = 39 
omp = 1
gpu = 0 
mem = 10000 
work_dir = './PE/' 
ff = GAFF2_mod() 

mol = utils.mol_from_smiles(smiles) 
ter = utils.mol_from_smiles(ter_smiles)
dp = poly.calc_n_from_num_atoms(mol, 1000, terminal1=ter)
homopoly = poly.polymerize_rw(mol, dp, tacticity='atactic') 
homopoly = poly.terminate_rw(homopoly, ter)
result = ff.ff_assign(homopoly, charge='gasteiger')
if not result: 
  print ('[ERROR: Can not assign force field parameters.]') 
# Generate simulation cell 
ac = poly.amorphous_cell(homopoly, 10, density=0.05)
# Equilibration MD 
eqmd = eq.EQ21step(ac, work_dir=work_dir) 
ac = eqmd.exec(temp=temp, press=1.0, mpi=mpi, omp=omp, gpu=gpu)
analy = eqmd.analyze() 
prop_data = analy.get_all_prop(temp=temp, press=1.0, save=True) 
result = analy.check_eq() 
# Additional equilibration MD 
for i in range(4): 
    if result: break 
    eqmd = eq.Additional(ac, work_dir=work_dir) 
    ac = eqmd.exec(temp=temp, press=press, mpi=mpi, omp=omp, gpu=gpu) 
    analy = eqmd.analyze() 
    prop_data = analy.get_all_prop(temp=temp, press=press, save=True) 
    result = analy.check_eq() 
if not result: 
  print('[ERROR: Did not reach an equilibrium state.]') 
# Non-equilibrium MD for thermal condultivity 
else: 
  nemd = tc.NEMD_MP(ac, work_dir=work_dir) 
  ac = nemd.exec(decomp=True, temp=temp, mpi=mpi, omp=omp, gpu=gpu) 
  nemd_analy = nemd.analyze() 
  TC = nemd_analy.calc_tc(decomp=True, save=True) 
  if not nemd_analy.Tgrad_data['Tgrad_check']: 
    print('[ERROR: Low linearity of temperature gradient.]') 
  print('Thermal conductivity: %f' % TC)
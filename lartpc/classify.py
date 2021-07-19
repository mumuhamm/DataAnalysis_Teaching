import ROOT
from larcv import larcv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 




chain_image2d = ROOT.TChain('image2d_data_tree')
chain_particle = ROOT.TChain('particle_mctruth_tree')
chain_image2d.AddFile('/eos/user/m/mumuhamm/test_40k.root')
chain_particle.AddFile('/eos/user/m/mumuhamm/test_40k.root')
print(chain_image2d.GetEntries(),'entries found!')

#get specific event


chain_image2d.GetEntry(0)
cpp_object = chain_image2d.image2d_data_branch
print('Object type:',cpp_object)

chain_particle.GetEntry(0)
cpp_object_particle = chain_particle.particle_mctruth_branch
print('Object type: {}\n'.format(str(cpp_object_particle)))

image2d_array = cpp_object.as_vector()


fig, axes = plt.subplots(1, image2d_array.size(), figsize=(12,4), facecolor='w')
for index,image2d in enumerate(image2d_array):
	image2d_numpy = larcv.as_ndarray(image2d)
	axes[index].imshow(image2d_numpy, interpolation='none',cmap='jet')
	nz_pixels=np.where(image2d_numpy>0.0)
	ylim = (np.min(nz_pixels[0])-5,np.max(nz_pixels[0])+5)
	xlim = (np.min(nz_pixels[1])-5,np.max(nz_pixels[1])+5)
	ylim = (np.max((ylim[0],0)), np.min((ylim[1],image2d_numpy.shape[1]-1)))
	xlim = (np.max((xlim[0],0)), np.min((xlim[1],image2d_numpy.shape[0]-1)))
	axes[index].set_ylim(ylim)
	axes[index].set_xlim(xlim)
plt.savefig("/eos/user/m/mumuhamm/www/mylar_40k.pdf")


print('Checking particle information for 1st entry...')
for particle in cpp_object_particle.as_vector():
	print('PDG Code: {}'.format(particle.pdg_code()))
	print('Initial energy: {:.3} GeV'.format(particle.energy_init()))




pdg_array      = np.zeros([chain_particle.GetEntries()],dtype=np.int32)
energy_array   = np.zeros([chain_particle.GetEntries()],dtype=np.float64)
momentum_array = np.zeros([chain_particle.GetEntries()],dtype=np.float64)


for index in range(chain_particle.GetEntries()):
	chain_particle.GetEntry(index)
	particle = chain_particle.particle_mctruth_branch.as_vector().front()
	pdg = int(particle.pdg_code())
	total_energy   = particle.energy_init() * 1000.
	kinetic_energy = total_energy - larcv.ParticleMass(pdg)
	momentum = np.sqrt(np.power(total_energy,2) - np.power(larcv.ParticleMass(pdg),2))

	pdg_array[index]      = pdg
	energy_array[index]   = kinetic_energy
	momentum_array[index] = momentum

	if momentum < 800:
		print(pdg,kinetic_energy,momentum)

df = pd.DataFrame(data={'pdg' : pdg_array, 'energy' : energy_array, 'momentum' : momentum_array})
pdg_list, pdg_counts = np.unique(df.pdg.values,return_counts=True)
print('PDGs found: {}'.format(pdg_list))
print('PDG counts: {}'.format(pdg_counts))

PDG2NAME = {11   : 'electron',
	22   : 'gamma',
	13   : 'muon',
	211  : 'pion',
	2212 : 'proton'}
for pdg in pdg_list:
	sub_df = df.query('pdg=={}'.format(pdg))
	min_value = sub_df.momentum.values.min()
	max_value = sub_df.momentum.values.max()
	print('{:10s} momentum range: {:.3g} => {:.3g} MeV/c'.format(PDG2NAME[pdg], min_value, max_value))


fig, ax = plt.subplots(figsize=(12,8),facecolor='w')

pdg_list = [13,22]
sub_df_E = df.query('pdg in {}'.format(pdg_list))
min_value = sub_df_E.energy.values.min()
max_value = sub_df_E.energy.values.max()


for pdg in pdg_list:
	pdg_df = sub_df_E.query('pdg == {}'.format(pdg))
	values = pdg_df.energy.values
	plt.hist(values, bins=40, range=(min_value,max_value), label='PDG {}'.format(pdg), alpha=0.5)

plt.tick_params(labelsize=20)
plt.grid(True,which='both')
plt.xlabel('Initial Kinetic Energy [MeV]',fontsize=20,fontweight='bold',fontname='monospace')
plt.ylabel('Number of entries',fontsize=20,fontweight='bold',fontname='monospace')
leg=plt.legend(fontsize=16,loc=4)
leg_frame=leg.get_frame()
leg_frame.set_facecolor('white')
plt.savefig("/eos/user/m/mumuhamm/www/photon_muon_energy.pdf")


fig, axm = plt.subplots(figsize=(12,8),facecolor='w')

pdg_list = [13,22]
sub_df_mom = df.query('pdg in {}'.format(pdg_list))
min_value_mom = sub_df_mom.momentum.values.min()
max_value_mom = sub_df_mom.momentum.values.max()

for pdg in pdg_list:
	pdg_df_mom = sub_df_mom.query('pdg == {}'.format(pdg))
	values = pdg_df_mom.momentum.values
	plt.hist(values, bins=40, range=(min_value_mom,max_value_mom), label='PDG {}'.format(pdg), alpha=0.5)

plt.tick_params(labelsize=20)
plt.grid(True,which='both')
plt.xlabel('Initial Momentum [MeV/c]',fontsize=20,fontweight='bold',fontname='monospace')
plt.ylabel('Number of entries',fontsize=20,fontweight='bold',fontname='monospace')
leg_mom=plt.legend(fontsize=16,loc=4)
leg_frame_mom=leg_mom.get_frame()
leg_frame_mom.set_facecolor('white')
plt.savefig("/eos/user/m/mumuhamm/www/photon_muon_momentum.pdf")


fig, axn = plt.subplots(figsize=(12,8),facecolor='w')

pdg_list = [2212,22]
sub_df_mom2 = df.query('pdg in {}'.format(pdg_list))
min_value_mom2 = sub_df_mom2.momentum.values.min()
max_value_mom2 = sub_df_mom2.momentum.values.max()


for pdg in pdg_list:
	pdg_df_mom2 = sub_df_mom2.query('pdg == {}'.format(pdg))
	values = pdg_df_mom2.momentum.values
	plt.hist(values, bins=40, range=(min_value_mom2,max_value_mom2), label='PDG {}'.format(pdg), alpha=0.5)

plt.tick_params(labelsize=20)
plt.grid(True,which='both')
plt.xlabel('Initial Momentum [MeV/c]',fontsize=20,fontweight='bold',fontname='monospace')
plt.ylabel('Number of entries',fontsize=20,fontweight='bold',fontname='monospace')
leg_mom2=plt.legend(fontsize=16,loc=2)
leg_frame_mom2=leg_mom2.get_frame()
leg_frame_mom2.set_facecolor('white')
plt.savefig("/eos/user/m/mumuhamm/www/photon_proton_momentum.pdf")

for index,image2d in enumerate(image2d_array):
	print(image2d.meta().dump())

sub_df_proton = df.query('pdg==2212 and energy < 50')
print('Found {} entries'.format(sub_df_proton.index.size))
print(sub_df_proton)

chain_image2d.GetEntry(32020)
cpp_object = chain_image2d.image2d_data_branch
print('Object type:',cpp_object)


image2dnew = cpp_object.as_vector().front()

fig, axl = plt.subplots(figsize=(12,12), facecolor='w')
image2d_numpy = larcv.as_ndarray(image2dnew)
axl.imshow(image2d_numpy, interpolation='none', vmin=0., vmax=image2d_numpy.mean(), cmap='gray')
# Find bounds for non-zero pixels + padding of 5 pixels
nz_pixels=np.where(image2d_numpy>0.0)
ylim = (np.min(nz_pixels[0])-5,np.max(nz_pixels[0])+5)
xlim = (np.min(nz_pixels[1])-5,np.max(nz_pixels[1])+5)
# Adjust for allowed image range
ylim = (np.max((ylim[0],0)), np.min((ylim[1],image2d_numpy.shape[1]-1)))
xlim = (np.max((xlim[0],0)), np.min((xlim[1],image2d_numpy.shape[0]-1)))
# Set range
axl.set_ylim(ylim)
axl.set_xlim(xlim)
plt.savefig("/eos/user/m/mumuhamm/www/pixel.pdf")


import ROOT
from ROOT import TChain
from larcv import larcv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
from matplotlib.ticker import FormatStrFormatter







def get_view_range(image2d):
	nz_pixels=np.where(image2d>0.0)
	ylim = (np.min(nz_pixels[0])-5,np.max(nz_pixels[0])+5)
	xlim = (np.min(nz_pixels[1])-5,np.max(nz_pixels[1])+5)
	ylim = (np.max((ylim[0],0)), np.min((ylim[1],image2d.shape[1]-1)))
	xlim = (np.max((xlim[0],0)), np.min((xlim[1],image2d.shape[0]-1)))
	return (xlim,ylim)

def show_event(entry):
	chain_image2d = ROOT.TChain('image2d_data_tree')
	chain_image2d.AddFile('test_10k.root')
	chain_label2d = ROOT.TChain('image2d_segment_tree')
	chain_label2d.AddFile('test_10k.root')
	print(chain_image2d.GetEntries(),'entries found!')
	if entry > 0:
		entry = np.random.randint(0,chain_label2d.GetEntries())
		print(entry)
		chain_label2d.GetEntry(entry)
		chain_image2d.GetEntry(entry)
		image2d = larcv.as_ndarray(chain_image2d.image2d_data_branch.as_vector().front())
		label2d = larcv.as_ndarray(chain_label2d.image2d_segment_branch.as_vector().front())
		xlim, ylim = get_view_range(image2d)
		fig, (ax0,ax1) = plt.subplots(1, 2, figsize=(10,10), facecolor='w')
		ax0.imshow(image2d, interpolation='none', cmap='jet', origin='lower')
		ax1.imshow(label2d, interpolation='none', cmap='jet', origin='lower',vmin=0., vmax=3.1)
		ax0.set_title('Data',fontsize=20,fontname='monospace',fontweight='bold')
		ax0.set_xlim(xlim)
		ax0.set_ylim(ylim)
		ax1.set_title('Label',fontsize=20,fontname='monospace',fontweight='bold')
		ax1.set_xlim(xlim)
		ax1.set_ylim(ylim)
		plt.savefig("/eos/user/m/mumuhamm/www/image_segment.pdf")
		return (np.array(image2d), np.array(label2d))

ENTRY = 10
image2d, label2d = show_event(ENTRY)

unique_values, unique_counts = np.unique(label2d, return_counts=True)
print('Label values:',unique_values)
print('Label counts:',unique_counts)

categories = ['Background','Shower','Track']
fig, axes = plt.subplots(1, len(unique_values), figsize=(10,6), facecolor='w')
xlim,ylim = get_view_range(image2d)

for index, value in enumerate(unique_values):
	ax = axes[index]
	mask = (label2d == value)
	ax.imshow(image2d * mask, interpolation='none', cmap='jet', origin='lower')
	ax.set_title(categories[index],fontsize=20,fontname='monospace',fontweight='bold')
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)
plt.savefig("/eos/user/m/mumuhamm/www/bkg_shower_track.pdf")


#for _ in range(5): show_event(ENTRY)


particle_mctruth_chain = ROOT.TChain("particle_mctruth_tree")
particle_mctruth_chain.AddFile("test_10k.root")
particle_mctruth_chain.GetEntry(ENTRY)
cpp_object = particle_mctruth_chain.particle_mctruth_branch

print('particle_mctruth_tree contents:')
for particle in cpp_object.as_vector():
	print(particle.dump())

particle_mcst_chain = TChain("particle_mcst_tree")
particle_mcst_chain.AddFile("test_10k.root")
particle_mcst_chain.GetEntry(ENTRY)
cpp_object_mcst = particle_mcst_chain.particle_mcst_branch

print('particle_mcst_tree contents:')
for particle in cpp_object_mcst.as_vector():
	print(particle.dump())



chain_image2dnew = ROOT.TChain('image2d_data_tree')
chain_image2dnew.AddFile('test_10k.root')
chain_image2dnew.GetEntry(ENTRY)

cpp_image2dnew = chain_image2dnew.image2d_data_branch.as_vector().front()
image2dnew = larcv.as_ndarray(cpp_image2dnew)

fig, axbox = plt.subplots(figsize=(8,6), facecolor='w')
axbox.imshow(image2dnew, interpolation='none', cmap='jet', origin='lower')
axbox.set_title('Data',fontsize=20,fontname='monospace',fontweight='bold')

for particle in cpp_object_mcst.as_vector():
	box = particle.boundingbox_2d().front()
	x_pixel_pos = box.min_x()  / cpp_image2dnew.meta().pixel_width()
	y_pixel_pos = box.min_y()  / cpp_image2dnew.meta().pixel_height()
	box_width   = box.width()  / cpp_image2dnew.meta().pixel_width()
	box_height  = box.height() / cpp_image2dnew.meta().pixel_height()
	axbox.add_patch(plt.Rectangle( (x_pixel_pos, y_pixel_pos), box_width, box_height,
		fill=False, 
		edgecolor='red', linewidth=3.5 )
		)
	axbox.text(x_pixel_pos, y_pixel_pos - 5,
			'PDG {:d} E={:.4g} MeV'.format(particle.pdg_code(),particle.energy_init()),
			bbox=dict(facecolor='yellow', alpha=0.5), fontsize=14, color='black'
			)


xlim, ylim = get_view_range(image2dnew)
axbox.set_xlim(xlim)
axbox.set_ylim(ylim)
plt.savefig("/eos/user/m/mumuhamm/www/particle_bounding_box.pdf")
		
num_entries = particle_mcst_chain.GetEntries()		
PDG_LIST = np.array([11,22,13,211,-211,2212])
PDG_NAME = ['electron','gamma','muon','piplus','piminus','proton']
total_deposit_energy  = np.zeros([num_entries],dtype=np.float32)
total_kinetic_energy  = np.zeros([num_entries],dtype=np.float32)
total_multiplicity    = np.zeros([num_entries],dtype=np.int32)
particle_multiplicity = [np.zeros([num_entries],dtype=np.int32) for _ in PDG_LIST]

for entry in np.arange(num_entries):
	particle_mcst_chain.GetEntry(entry)
	mcstobject = particle_mcst_chain.particle_mcst_branch
	deposit_energy = 0.
	initial_energy = 0.
	multiplicity = 0
	for particle in mcstobject.as_vector():
		deposit_energy += particle.energy_deposit()
		if not particle.track_id() == particle.parent_track_id(): continue
		particle_multiplicity[np.where(particle.pdg_code() == PDG_LIST)[0][0]][entry] += 1
		multiplicity += 1
		initial_energy += (particle.energy_init() - larcv.ParticleMass(particle.pdg_code()))
	total_multiplicity[entry]   = multiplicity
	total_deposit_energy[entry] = deposit_energy
	total_kinetic_energy[entry] = initial_energy

df = pd.DataFrame(data={'initial_energy' : total_kinetic_energy,
	 'deposit_energy' : total_deposit_energy,
	 'multi_total'    : total_multiplicity,
	 'multi_electron' : particle_multiplicity[0],
	 'multi_gamma'    : particle_multiplicity[1],
	 'multi_muon'     : particle_multiplicity[2],
	 'multi_piplus'   : particle_multiplicity[3],
	  'multi_piminus'  : particle_multiplicity[4],
	  'multi_proton'   : particle_multiplicity[5]})

unique_values, unique_counts = np.unique(df.multi_total, return_counts=True)
print()
print('Multiplicity values:',unique_values)
print('Multiplicity counts:',unique_counts)

fig, axes = plt.subplots(1,3,figsize=(16,4),facecolor='w')
for index, name in enumerate(['electron','gamma','muon']):
	ax = axes[index]
	exec('ax.hist(df.multi_{}, bins=5, range=(-0.5,4.5))'.format(name))
	ax.set_xlim(-0.5,4.5)
	ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
	ax.xaxis.set_ticks(np.arange(0, 5, 1.))
	ax.set_title(name,fontsize=20,fontweight='bold',fontname='monospace')
	ax.set_xlabel('Multiplicity',fontsize=20,fontweight='bold',fontname='monospace')
	ax.tick_params('x',labelsize=20)
	ax.tick_params('y',labelsize=16)
	ax.grid()
plt.savefig("/eos/user/m/mumuhamm/www/multiplicity_e_gam_muon.pdf")
print('Event counts with at least 1 electron or muon:',df.query('multi_electron >0 or multi_muon >0').index.size)
print('Event counts with at least 1 gamma ray:',df.query('multi_gamma>0').index.size)



fig, ax = plt.subplots(figsize=(8,6),facecolor='w')

track_counts  = df.multi_muon.values + df.multi_piplus.values + df.multi_piminus.values + df.multi_proton.values
shower_counts = df.multi_electron.values + df.multi_gamma.values

ax.hist(track_counts,  bins=6, range=(-0.5,5.5), alpha=0.5, label='track multiplicity')
ax.hist(shower_counts, bins=6, range=(-0.5,5.5), alpha=0.5, label='shower multiplicity')

ax.set_xlim(-0.5,5.5)
ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
ax.xaxis.set_ticks(np.arange(0, 6, 1.))
ax.set_title("Track vs. Shower Multiplicity",fontsize=20,fontweight='bold',fontname='monospace')
ax.set_xlabel('Multiplicity',fontsize=20,fontweight='bold',fontname='monospace')
ax.tick_params('x',labelsize=20)
ax.tick_params('y',labelsize=16)
ax.grid()

leg=plt.legend(fontsize=16,loc=1)
leg_frame=leg.get_frame()
leg_frame.set_facecolor('white')
plt.savefig("/eos/user/m/mumuhamm/www/trackvsshower_multiplicity.pdf")



fig,ax = plt.subplots(figsize=(10,9),facecolor='w')
plt.hist(df.initial_energy.values,bins=40,range=(0.,4000),alpha=0.5,label="Initial Kinetic Energy")
plt.hist(df.deposit_energy.values,bins=40,range=(0.,4000),alpha=0.5,label="Deposited Energy")
plt.tick_params(labelsize=20)
plt.grid(True,which='both')
plt.xlabel('Energy [MeV]',fontsize=20,fontweight='bold',fontname='monospace')
plt.ylabel('Number of entries',fontsize=20,fontweight='bold',fontname='monospace')
leg=plt.legend(fontsize=20,loc=1)
leg_frame=leg.get_frame()
leg_frame.set_facecolor('white')
plt.savefig("/eos/user/m/mumuhamm/www/energydistribution.pdf")


fig, axes = plt.subplots(1,3,figsize=(16,6),facecolor='w')

queries = ['multi_total == (multi_electron + multi_gamma)',
		 'multi_total == (multi_proton)',
		  'multi_total == (multi_muon + multi_piplus + multi_piminus + multi_proton)']
titles  = ['Electron + Gamma','Proton','Muon + Pion']

for index, query in enumerate(queries):
	sub_df=df.query(query)
	ax = axes[index]
	ax.hist(sub_df.initial_energy.values,bins=25,range=(0.,4000),alpha=0.5,label='Initial')
	ax.hist(sub_df.deposit_energy.values,bins=25,range=(0.,4000),alpha=0.5,label='Deposited')
	ax.set_title(titles[index],fontsize=14,fontweight='bold',fontname='monospace')
	ax.tick_params(labelsize=14)
	ax.grid(True,which='both')
	ax.set_xlabel('Energy [MeV]',fontsize=14,fontweight='bold',fontname='monospace')
	ax.set_ylabel('Number of entries',fontsize=14,fontweight='bold',fontname='monospace')
	ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
	ax.xaxis.set_ticks(np.arange(0, 4000, 1000.))
	leg=ax.legend(fontsize=14,loc=1)
	leg_frame=leg.get_frame()
	leg_frame.set_facecolor('white')

plt.savefig("/eos/user/m/mumuhamm/www/three_categories.pdf")


chain_particle  = ROOT.TChain("particle_mcst_tree")
chain_cluster2d = ROOT.TChain("cluster2d_mcst_tree")

for chain in [chain_particle, chain_cluster2d]:
	chain.AddFile('test_10k.root')
	chain.GetEntry(ENTRY)

cpp_particle  = chain_particle.particle_mcst_branch
cpp_cluster2d = chain_cluster2d.cluster2d_mcst_branch


for index, particle in enumerate(cpp_particle.as_vector()):
	cluster = cpp_cluster2d.cluster_pixel_2d(0).as_vector()[index]
	msg = 'Particle {:d} has TrackID {:d} and PDG {:d} with {:d} non-zero pixels!'
	print(msg.format(index,particle.track_id(),particle.pdg_code(),cluster.as_vector().size()))

num_particles = cpp_particle.as_vector().size()
fig, axes = plt.subplots(1,num_particles,figsize=(24,8), facecolor='w')
for index, particle in enumerate(cpp_particle.as_vector()):
	ax = axes[index]
	title = 'TrackID {:d} PDG {:d}'.format(particle.track_id(),particle.pdg_code())
	cluster2d = cpp_cluster2d.cluster_pixel_2d(0).as_vector()[index]
	image2d = larcv.as_image2d(cluster2d, cpp_cluster2d.cluster_pixel_2d(0).meta())
	image2d = larcv.as_ndarray(image2d) * 100.
	ax.imshow(image2d, interpolation='none', cmap='jet', origin='lower')
	ax.set_title(title,fontsize=14,fontname='monospace',fontweight='bold')
plt.savefig("/eos/user/m/mumuhamm/www/clusters.pdf")

cpp_image2d = chain_image2dnew.image2d_data_branch.as_vector().front()
fig, (ax0,ax1) = plt.subplots(1,2,figsize=(16,8), facecolor='w')
image2d = larcv.as_ndarray(cpp_image2d)
ax0.imshow(image2d, interpolation='none', cmap='jet', vmin=0, vmax=1000, origin='lower')
ax0.set_title('image2d_data_tree',fontsize=24,fontname='monospace',fontweight='bold')
cluster2d = np.zeros(image2d.shape,dtype=np.float32)
num_particles = cpp_particle.as_vector().size()
for index, particle in enumerate(cpp_particle.as_vector()):
	title = 'TrackID {:d} PDG {:d}'.format(particle.track_id(),particle.pdg_code())
	particle_cluster = cpp_cluster2d.cluster_pixel_2d(0).as_vector()[index]
	image2d_from_cluster2d = larcv.as_image2d(particle_cluster, cpp_cluster2d.cluster_pixel_2d(0).meta())
	image2d_from_cluster2d = larcv.as_ndarray(image2d_from_cluster2d) * 100.
	cluster2d = cluster2d + image2d_from_cluster2d

ax1.imshow(cluster2d, interpolation='none', cmap='jet', vmin=0, vmax=1000, origin='lower')
ax1.set_title('cluster2d_mcst_tree',fontsize=24,fontname='monospace',fontweight='bold')
plt.savefig("/eos/user/m/mumuhamm/www/clusters_and_2dimage_comparison.pdf")


msg = 'Image from {:s} ... mean {:.4g} pixel count {:d}'
print(msg.format('image2d_data_tree',image2d.mean(),(image2d>0).astype(np.int32).sum()))
print(msg.format('cluster2d_data_tree',cluster2d.mean(),(cluster2d>0).astype(np.int32).sum()))


msg = 'Image from {:s} ... mean {:.4g} pixel count {:d}'
nonzero_image2d   = image2d[np.where(image2d>=10.)]
nonzero_cluster2d = cluster2d[np.where(cluster2d>=10.)]
print(msg.format('image2d_data_tree',nonzero_image2d.mean(),len(nonzero_image2d)))
print(msg.format('cluster2d_data_tree',nonzero_cluster2d.mean(),len(nonzero_cluster2d)))


fig, (ax0,ax1) = plt.subplots(1,2,figsize=(16,8), facecolor='w')
ax0.imshow(cluster2d, interpolation='none', cmap='jet', vmin=0, vmax=1000, origin='lower')
ax0.set_title('Trajectory',fontsize=24,fontname='Georgia',fontweight='bold')

label_cluster2d = np.zeros(cluster2d.shape,dtype=np.float32)
num_particles = cpp_particle.as_vector().size()

for index, particle in enumerate(cpp_particle.as_vector()):
	title = 'TrackID {:d} PDG {:d}'.format(particle.track_id(),particle.pdg_code())
	particle_cluster = cpp_cluster2d.cluster_pixel_2d(0).as_vector()[index]
	image2d_from_cluster2d = larcv.as_image2d(particle_cluster, cpp_cluster2d.cluster_pixel_2d(0).meta())
	image2d_from_cluster2d = larcv.as_ndarray(image2d_from_cluster2d) * 100.
	label_cluster2d[np.where(image2d_from_cluster2d >= 10.)] = index+1

ax1.imshow(label_cluster2d, interpolation='none', cmap='jet', vmin=0, vmax=num_particles, origin='lower')
ax1.set_title('Instance labels',fontsize=24,fontname='monospace',fontweight='bold')

plt.savefig("/eos/user/m/mumuhamm/www/instance_wise_level.pdf")





































































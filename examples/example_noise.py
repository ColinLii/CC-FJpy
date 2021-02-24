import matplotlib.pyplot as plt
import numpy as np
import ccfj


data = np.load('summed.npz')
ncfs = data['ncfs']
nc = 1000
c = np.linspace(2000,5000,nc)
r = data['r']*1e3
f = data['f']

ds00 = ccfj.fj_noise(np.real(ncfs),r,c,f,fstride=1,itype=0,func=0)
ds01 = ccfj.fj_noise(np.real(ncfs),r,c,f,fstride=1,itype=1,func=0)
ds10 = ccfj.fj_noise(np.real(ncfs),r,c,f,fstride=1,itype=0,func=1)
ds11 = ccfj.fj_noise(np.real(ncfs),r,c,f,fstride=1,itype=1,func=1)

fig,ax=plt.subplots(nrows=2,ncols=2,figsize=(10,10))
ax[0][0].pcolormesh(f,c/1e3,ds00,cmap='jet',vmin=0,vmax=0.8)
ax[0][0].set_xlim([0,0.5])
ax[0][0].set_ylabel('Phase velocity (km/s)')
ax[0][0].set_xlabel('Frequency (Hz)')

ax[0][1].pcolormesh(f,c/1e3,ds01,cmap='jet',vmin=0,vmax=0.8)
ax[0][1].set_xlim([0,0.5])
ax[0][1].set_ylabel('Phase velocity (km/s)')
ax[0][1].set_xlabel('Frequency (Hz)')

ax[1][0].pcolormesh(f,c/1e3,ds10,cmap='jet',vmin=0,vmax=0.8)
ax[1][0].set_xlim([0,0.5])
ax[1][0].set_ylabel('Phase velocity (km/s)')
ax[1][0].set_xlabel('Frequency (Hz)')

ax[1][1].pcolormesh(f,c/1e3,ds11,cmap='jet',vmin=0,vmax=0.8)
ax[1][1].set_xlim([0,0.5])
ax[1][1].set_ylabel('Phase velocity (km/s)')
ax[1][1].set_xlabel('Frequency (Hz)')

plt.show()

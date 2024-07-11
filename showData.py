import numpy as np
import matplotlib.pyplot as plt


baVector = np.linspace(0, 6, 40)
getData = np.zeros((40, 2))
temp = 0

for l in baVector:
	getData[temp, 0] = l

	cmd = '/home/weicheng/Desktop/seu_project/magnetic/3D/JMPS_61/datafiles/simDER_rodLength_0.250000_Ba_' + '%7.6f' % l + '.txt'
	data123 = np.loadtxt(cmd)
	getData[temp, 1] = data123


	#cmd = '/home/weicheng/Desktop/seu_project/magnetic/3D/JMPS_61/datafiles/simDER_rodLength_0.125000_Ba_' + '%7.6f' % l + '.txt'
	#data123 = np.loadtxt(cmd)
	#getData[temp, 2] = data123

	#cmd = '/home/weicheng/Desktop/seu_project/magnetic/3D/JMPS_61/datafiles/simDER_rodLength_0.150000_Ba_' + '%7.6f' % l + '.txt'
	#data123 = np.loadtxt(cmd)
	#getData[temp, 2] = data123


	#cmd = '/home/weicheng/Desktop/seu_project/magnetic/3D/JMPS_61/datafiles/simDER_rodLength_0.175000_Ba_' + '%7.6f' % l + '.txt'
	#data123 = np.loadtxt(cmd)
	#getData[temp, 4] = data123
	
	#cmd = '/home/weicheng/Desktop/seu_project/magnetic/3D/JMPS_61/datafiles/simDER_rodLength_0.200000_Ba_' + '%7.6f' % l + '.txt'
	#data123 = np.loadtxt(cmd)
	#getData[temp, 3] = data123


	temp = temp + 1

#getData[:,1] = np.absolute(getData[:,1] - getData[0, 1])
#getData[:,2] = np.absolute(getData[:,2] - getData[0, 2])
#getData[:,3] = np.absolute(getData[:,3] - getData[0, 3])

plt.plot(getData[:,0], getData[:,1], '-o')

np.savetxt('/home/weicheng/Desktop/seu_project/magnetic/3D/JMPS_61/saveData/l_0.25.txt', getData)

#plt.plot(getData[:,0], getData[:,2], '-s')
#plt.plot(getData[:,0], getData[:,3], '-^')

#np.savetxt('/home/weicheng/Desktop/seu_project/magnetic/3D/JMPS_61/saveData/buckling_l.txt',getData)

'''
getData = np.zeros((400, 3))

baVector = np.logspace(-1, 1, 20)
RodLength = np.linspace(0.1, 0.5, 20)

temp = 0
for l in baVector:
	for r in RodLength:
		cmd = '/home/weicheng/Desktop/seu_project/magnetic/3D/JMPS_61/datafiles_2/simDER_rodLength_' + '%7.6f' % r + '_Ba_' + '%7.6f' % l + '.txt'
		data123 = np.loadtxt(cmd)

		getData[temp, 0] = r
		getData[temp, 1] = l

		if data123 > 0.5:
			getData[temp, 2] = 1
			plt.plot(r, l, 'rs')

		if data123 < 0.5:
			getData[temp, 2] = 0
			plt.plot(r, l, 'bo')

		temp = temp + 1

x = np.linspace(0.1, 0.5 ,100)
y = np.linspace(0.0, 0.0 ,100)

temp = 0
for aa in x:
	y[temp] = 0.0423 / (aa * aa)
	temp = temp + 1

plt.plot(x, y, '-')

plt.yscale('log')

np.savetxt('/home/weicheng/Desktop/seu_project/magnetic/3D/JMPS_61/saveData/data.txt', getData)
'''

plt.show()
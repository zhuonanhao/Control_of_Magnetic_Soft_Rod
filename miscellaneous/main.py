import numpy as np

f = open('Commands.txt','w')

#baVector = np.logspace(-1, 1, 20)

baVector = np.linspace(0, 6, 40)
RodLength = np.linspace(0.25, 0.25, 1)

print(baVector)

for l in baVector:
	for r in RodLength:
		cmdline = 'export OMP_NUM_THREADS=1; ./simDER' + ' option.txt ' + '--' + ' RodLength ' + '%7.6f' % r + ' baVector ' + '%7.6f' % l + ' 0.0 0.0' + '\n'
		f.write(cmdline)
	
'''
for l in baVector:
	cmdline = 'export OMP_NUM_THREADS=1; ./simDER' + ' option.txt ' + '--' + ' RodLength 0.15 baVector ' + '%7.6f' % l + ' 0.0 0.0' + '\n'
	f.write(cmdline)

for l in baVector:
	cmdline = 'export OMP_NUM_THREADS=1; ./simDER' + ' option.txt ' + '--' + ' RodLength 0.2 baVector ' + '%7.6f' % l + ' 0.0 0.0' + '\n'
	f.write(cmdline)

'''

f.close()

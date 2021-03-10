import numpy as np
import matplotlib.pyplot as plt

def knee_center():

	c1 = 0.0791
	c2 = 5.733E-4
	c3 = 7.682E-6
	c4 = 5.759E-8

	c5 = 0.3695
	c6 = 2.958E-3
	c7 = 7.666E-6

	c8 = -0.0683
	c9 = 8.804E-4
	c10 = 3.750E-6

	c11 = -0.1283
	c12 = 4.796E-4

	for Beta in range(0, 120):

		Varus = (c1 * Beta) - (c2 * Beta**2) - (c3 * Beta**3) + (c4 * Beta**4)
		IntRot = (c5 * Beta) - (c6 * Beta**2) + (c7 * Beta**3)
		ydist = (c8 * Beta) + (c9 * Beta**2) - (c10 * Beta**3)
		zdist = (c11 * Beta) + (c12 * Beta**2)

		X1 = -60
		Y1 = -np.sin(Varus) * X1 + ydist
		Z1 = np.cos(Varus) * np.sin(IntRot) * X1 + zdist

def _FBM_angles(L1, L2, L3, L4, theta1, theta3):

	K = (L1**2 + L2**2 + L3**2 + L4**2)/2
	A = K - L2*L3*np.cos(theta3) + L1*L2*np.cos(theta1) + L1*L3*np.cos(theta1-theta3)
	B = 2*(L2*L3*np.sin(theta3) - L1*L2*np.sin(theta1))
	C = K + L2*L3*np.cos(theta3) - L1*L2*np.cos(theta1) - L1*L3*np.cos(theta1-theta3)
	
	if B**2 + 4*A*C < 0:
		config["draw":0]
		return config
	
	theta21 = 2*np.arctan((-B+(B**2 + 4*A*C)**0.5)/(2*A))
	theta22 = 2*np.arctan((-B-(B**2 + 4*A*C)**0.5)/(2*A))
	
	
	if -1<=((L2*np.sin(theta21)+ L3*np.sin(theta3) - L1*np.sin(theta1))/L4)<=1:
		theta4 = np.arcsin((L2*np.sin(theta21)+ L3*np.sin(theta3) - L1*np.sin(theta1))/L4)
		theta2 = theta21
	else:
		theta4 = np.arcsin((L2*np.sin(theta22)+ L3*np.sin(theta3) - L1*np.sin(theta1))/L4)
		theta2 = theta22

	config = {'L1':L1,'L2':L2,'L3':L3,'L4':L4,
			'theta1':theta1,'theta2':theta2,
			'theta3':theta3,'theta4':theta4,
			'draw':1}
	return config
	
def FBM_simulations(config):

	config['X7'] =  config['X5'] + config['L3']*np.cos(-1*config['theta3'])
	config['X8'] =  config['X6'] + config['L3']*np.sin(-1*config['theta3'])

	config['X1'] =  config['X7'] - config['L4']*np.cos(config['theta4'])
	config['X2'] =  config['X8'] - config['L4']*np.sin(config['theta4'])

	config['X3'] =  config['X1'] - config['L1']*np.cos(-1*config['theta1'])
	config['X4'] =  config['X2'] - config['L1']*np.sin(-1*config['theta1'])


def plot(config):
	print(config)
	l1=ax.plot([config['X1'], config['X3']],[config['X2'], config['X4']])
	l2=ax.plot([config['X3'], config['X5']],[config['X4'], config['X6']])
	l3=ax.plot([config['X5'], config['X7']],[config['X6'], config['X8']])
	l4=ax.plot([config['X7'], config['X1']],[config['X8'], config['X2']])
	plt.pause(0.1)
	l1.pop(0).remove()
	l2.pop(0).remove()
	l3.pop(0).remove()
	l4.pop(0).remove()
	
fig,ax = plt.subplots(1,1, figsize=(5,5))
for i in range(36):
	config = _FBM_angles(1,1,1,1,np.pi/6, 3*np.pi/4 + i*np.pi/18)
	if not config['draw']:
		print("Error")
		continue
	config['X5'] = 0
	config['X6'] = 0
	FBM_simulations(config)
	plot(config)

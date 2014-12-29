"""
Author: Jian Dai
"""
#!/usr/bin/python
from math import sin, cos, pi, acos
import numpy as np

class CalcInterVecAng:
	def __init__(self, infilename, outfilename):
		self.outfilename = outfilename
		phipsi_lst = []
		for line in open(infilename, 'r'):
			resid1 = int(line.split()[0])
			phi = float(line.split()[2])
			psi = float(line.split()[3])
			phipsi_lst.append([phi, psi, resid1])
		self.phipsi_lst = phipsi_lst

	def R1(self, angle):
		angle = angle * pi/180.0
		ret_matrix = np.array([[1.0, 0.0, 0.0], \
			                   [0.0, cos(angle), -1.0*sin(angle)], \
			                   [0.0, sin(angle), cos(angle)]])
		return ret_matrix

	def R3(self, angle):
		angle = angle * pi/180.0
		ret_matrix = np.array([[cos(angle), -1.0*sin(angle), 0.0], \
			                   [sin(angle), cos(angle), 0.0], \
			                   [0.0, 0.0, 1.0]])
		return ret_matrix

	def calcOneInterVecAng(self, i, j):
		resid_i = self.phipsi_lst[i][2]
		resid_j = self.phipsi_lst[j][2]
		# transformation matrix from DFF to PAF
		tt = np.array([[0.719, 0.0, -0.695], \
			           [-0.695, 0.0, -0.719], \
			           [0.0, 1.0, 0.0]])
		extAlpha = 65.0
		extBeta = 59.0
		extGamma = 70.0
		extOmega = 180.0
		betaNHRad = 17.0 * pi/180.0
		R3gamma = self.R3(extGamma)
		R3alpha = self.R3(extAlpha)
		R1omega = self.R1(extOmega)
		R3beta = self.R3(extBeta)
		# Coorindate of NH vector in the DFF
		DFF_NH0 = np.array([sin(betaNHRad), 0.0, cos(betaNHRad)])
		PAF_NH0 = tt.dot(DFF_NH0)

		R = np.eye(3)
		for k in range(i, j):
			phi = self.phipsi_lst[k][0]
			psi = self.phipsi_lst[k][1]
			R1phi = self.R1(phi)
			R1psi = self.R1(psi)
			R = R.dot(R1phi).dot(R3gamma).dot(R1psi).dot(R3alpha).dot(R1omega).dot(R3beta)

		PAF_NH = R.dot(PAF_NH0)
		dotprod = PAF_NH0.dot(PAF_NH)
		interVecAng = acos(dotprod) * 180.0/pi
		return [resid_i, resid_j, interVecAng]

	def calcAllInterVecAng(self):
		outfile = open(self.outfilename, 'w')
		N = len(self.phipsi_lst)
		counter = 0
		for i in range(0, N-1):
			for j in range(i+1, N):
				[resid_i, resid_j, interVecAng] = self.calcOneInterVecAng(i, j)
				print>>outfile, "%4i%4i%4i%8.3f" % (counter, resid_i, resid_j, interVecAng)
				counter += 1
		outfile.close()

if __name__ == '__main__':
	#instance = CalcInterVecAng('ideal_phipsi.dat', 'ideal_intervec_ang_pisema.xvg')
	instance = CalcInterVecAng('phipsi.dat', 'intervec_ang.xvg')
	instance.calcAllInterVecAng()

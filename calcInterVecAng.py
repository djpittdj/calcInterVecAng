#!/usr/bin/python
"""
Author: Jian Dai
"""
from math import sin, cos, pi, acos
import numpy as np

class CalcInterVecAng:
	def __init__(self, infilename, outfilename, isStruct=False):
		self.outfilename = outfilename
		phipsi_lst = []
		for line in open(infilename, 'r'):
			resid1 = int(line.split()[0])
			phi = float(line.split()[2])
			psi = float(line.split()[3])
			phipsi_lst.append([phi, psi, resid1])
		self.phipsi_lst = phipsi_lst
		self.isStruct = isStruct

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
		# transformation matrix from PAF to DFF, page 53 of Jun Hu's thesis
		# PAF = DFF.M1 -> DFF = PAF.M1_inv

		if self.isStruct: # Calculate from calcExtAngle.tcl on Chimera structure
			M1 = self.R3(-43.2).dot(self.R1(90.0))
			extAlpha = 63.4
			extBeta = 58.2
			extGamma = 70.0
			betaNHRad = 14.0 * pi/180.0
		else: # Using ideal peptide plane geometry
			M1 = self.R3(-44).dot(self.R1(90.0))
			extAlpha = 65.0
			extBeta = 59.0
			extGamma = 70.0
			betaNHRad = 17.0 * pi/180.0

		M1_inv = np.linalg.pinv(M1)
		extOmega = 180.0

		R3gamma = self.R3(extGamma)
		R3alpha = self.R3(extAlpha)
		R1omega = self.R1(extOmega)
		R3beta = self.R3(extBeta)
		# Coorindate of NH vector in the PAF
		PAF_NH0 = np.array([sin(betaNHRad), 0.0, cos(betaNHRad)])
		DFF_NH0 = PAF_NH0.dot(M1_inv)

		R = np.eye(3)
		for k in range(i, j):
			phi = self.phipsi_lst[k][0]
			psi = self.phipsi_lst[k][1]
			R1phi = self.R1(phi)
			R1psi = self.R1(psi)
			# page 54 of Jun Hu's thesis
			R = R.dot(R1phi).dot(R3gamma).dot(R1psi).dot(R3alpha).dot(R1omega).dot(R3beta)

		DFF_NH = DFF_NH0.dot(R)
		dotprod = DFF_NH0.dot(DFF_NH)
		interVecAng = acos(dotprod) * 180.0/pi
		return (resid_i, resid_j, interVecAng)

	def calcAllInterVecAng(self):
		outfile = open(self.outfilename, 'w')
		N = len(self.phipsi_lst)
		counter = 0
		for i in range(0, N-1):
			for j in range(i+1, N):
				(resid_i, resid_j, interVecAng) = self.calcOneInterVecAng(i, j)
				print>>outfile, '%4i%4i%4i%8.3f' % (counter, resid_i, resid_j, interVecAng)
				counter += 1
		outfile.close()

if __name__ == '__main__':
	instance = CalcInterVecAng('ideal_phipsi.dat', 'ideal_intervec_ang_pisema.xvg', isStruct=True)
	instance.calcAllInterVecAng()
	instance = CalcInterVecAng('phipsi.dat', 'intervec_ang.xvg', isStruct=False)
	instance.calcAllInterVecAng()

import numpy as np


class Molecule:

    def __init__(self, xyzfile):

        self.xyzfile = xyzfile
        self.atoms = []
        self.coords = []
        with open(self.xyzfile, 'r') as fp:
            self.natom = int(fp.readline())
            self.title = fp.readline().strip()
            for line in fp:
                atom, x, y, z = line.split()
                self.atoms.append(int(atom))
                self.coords.append([float(x), float(y), float(z)])

        self.coords = np.array(self.coords)
        self.occ = sum(self.atoms)/2
        

    def print_geom(self):

        if self.title:
            print(f'Name of the geometry: {self.title}')
        print(f'Number of Atoms : {self.natom}')
        print('Coordinate File: ')
        for i, j in zip(self.atoms, self.coords):
            print(f'{i}{j[0]:12.6f}{j[1]:12.6f}{j[2]:12.6f}')


    def bond(self, i, j):

        return np.linalg.norm(self.coords[i]- self.coords[j])


    def bond_angle(self, i, j, k):

        ji = self.coords[i] - self.coords[j]
        jk = self.coords[k] - self.coords[j] 

        cosine_angle = np.dot(ji, jk)/(np.linalg.norm(ji)*np.linalg.norm(jk))
        angle = np.arccos(cosine_angle)
        return np.degrees(angle)


    def dihedral_angle(self, a, b, c, d):
        """
        Calculates the dihedral angle between four points a, b, c, and d.
        """
        a = self.coords[a]
        b = self.coords[b]
        c = self.coords[c]
        d = self.coords[d]

        b1 = b - a
        b2 = c - b
        b3 = d - c

        # normalize b2
        b2 /= np.linalg.norm(b2)

        # calculate the normal vectors of the planes defined by the three successive points
        v1 = np.cross(b1, b2)
        v2 = np.cross(b2, b3)

        # calculate the dihedral angle between the planes
        angle = np.arctan2(np.linalg.norm(np.cross(v1, v2)), np.dot(v1, v2))

        return np.degrees(angle) 


    # def com(self):
    #     masses = {6:12.0000000, 8:15.994914619, 1:1.00782503223}

    #     # Calculate Molecular Mass
    #     M = 0
    #     for i in range(len(self.atoms)):
    #         M += masses[int(self.atoms[i])]

    #     print(f'Molecular Mass of the Molecule is {M:6f}')


    #     xcm = 0
    #     ycm = 0
    #     zcm = 0

    #     for i, j in zip(self.atoms, self.coords):
    #         # print(i, masses[int(i)]) 
    #         xcm += masses[int(i)]*j[0]
    #         ycm += masses[int(i)]*j[1]
    #         zcm += masses[int(i)]*j[2]

    #     xcm = xcm/M
    #     ycm = ycm/M
    #     zcm = zcm/M 

    #     return xcm, ycm, zcm


    def translate(self, x, y, z):

        coords = self.coords
        for i in range(len(self.coords)):
            coords[i][0] += x
            coords[i][1] += y
            coords[i][2] += z
        return coords

    def read_hessian(self, hess_file):


        # Reading Manually
        # with open(hess_file, 'r') as fp:
        #     natoms = int(fp.readline().strip())
        #     H = [[0.0]*(natoms*3) for i in range(natoms*3)]
        #     for i in range(natoms*3):
        #         for j in range(natoms):
        #             row = fp.readline().split()
        #             H[i][3*j] = float(row[0])
        #             H[i][3*j+1] = float(row[1])
        #             H[i][3*j+2] = float(row[2])
        # H = np.array(H)

        # Reading with numpy

        with open(hess_file, 'r') as fp:
            H = np.array([line.split() for line in fp.readlines()[1:]], dtype=float).reshape(self.natom*3, self.natom*3)
        return H


    def read_intregals(self, intregal):
        with open(intregal, 'r') as fp:
            lines = fp.readlines()
        nao = int(lines[-1].strip().split()[0])
        matrix = np.zeros([nao, nao])
        for line in lines:
            i, j, val = line.strip().split()
            i, j, val = int(i) - 1, int(j) - 1, float(val)
            matrix[i][j] = matrix[j][i]  = val

        return matrix

    def read_2e_intregals(self, intregal):
        with open(intregal, 'r') as fp:
            lines = fp.readlines()
        nao = int(lines[-1].strip().split()[0])
        matrix = np.zeros([nao, nao, nao, nao])
        for line in lines:
            i, j, k, l, val = line.strip().split()
            i, j, k, l, val = int(i) - 1, int(j) - 1, int(k) - 1, int(l) - 1, float(val)
            matrix[i][j][k][l] = matrix[i][j][l][k] = matrix[j][i][k][l] = matrix[j][i][l][k] \
                = matrix[k][l][i][j] = matrix[k][l][j][i] = matrix[l][k][i][j] = matrix[l][k][j][i] = val 

        return matrix












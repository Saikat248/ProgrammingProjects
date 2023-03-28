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
                self.atoms.append(atom)
                self.coords.append([float(x), float(y), float(z)])
        self.coords = np.array(self.coords)
        

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





# mol = Molecule('geom.dat')

# print(mol.bond_angle(0,1,2))








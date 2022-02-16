# Processing functions library (pflibs)
import re

import numpy as np
from scipy.spatial.transform import Rotation


def rot_wcha_modified(data, axis, angle):
    for n in range(len(axis)):
        r = Rotation.from_euler(axis[n], angle[n], degrees=True).as_matrix()
        data = data @ r.T
    return data


def laue_transform(data, phi, chi, theta):
    Mr = rot_wcha_modified(data, ['y', 'z', 'x'], [-theta, chi - 90, -phi])
    in_plane = Mr @ np.array([1, 0, 0])
    out_plane = Mr @ np.array([0, 1, 0])

    return in_plane, out_plane


def get_specangles(file):
    phi = chi = th = gam = delta = None
    f = open(file)
    lines = f.readlines()
    for line in lines:
        if "phi" in line:
            phi = float(line.split("=")[1])
        if "chi" in line:
            chi = float(line.split("=")[1])
        if "th" in line:
            th = float(line.split("=")[1])
        if "gam" in line:
            gam = float(line.split("=")[1])
        if "del" in line:
            delta = float(line.split("=")[1])

    return th, chi, phi, delta, gam


# extract error
def extract_err():
    find_goodness = re.compile(r"[$]goodness\d\t+(\d+.*\d+)\t")
    find_rms = re.compile(r"rms_error\d\s+(.+?)\t")

    file_name = open("./Index.txt").read()
    file_split = re.split("[$]pattern\d", file_name)[1:]
    good_ness = []
    rms_list = []

    for patt in file_split:
        goodness = np.round(float(find_goodness.findall(patt)[0]), 4)
        rms_err = np.round(float(find_rms.findall(patt)[0]), 5)
        good_ness.append(goodness)
        rms_list.append(rms_err)
        print(goodness)
        print(rms_list)

    return good_ness, rms_list


def extract_rmat():
    find_rmat = re.compile(r"[$]rotation_matrix\d\t+{(.+)}")
    file_name = open("./Index.txt").read()
    file_split = re.split("[$]pattern\d", file_name)[1:]

    rmat_list = []
    for patt in file_split:
        # print(find_rmat.findall(patt))
        rmat_string = find_rmat.findall(patt)[0]
        # print(rmat_string)
        rmat_colum_str = re.findall(r"{[-+]?\d+\.?\d*,[-+]?\d+\.?\d*,[-+]?\d+\.?\d*}", rmat_string)
        # print(rmat_colum_str)
        rmat_colum_list = []
        for rmat_colum in rmat_colum_str:
            rmat_colum_list.append([float(s) for s in re.findall(r"[-+]?\d+\.?\d*", rmat_colum)])
        rmat_colum_array = np.array(rmat_colum_list)
        rmat_list.append(rmat_colum_array)

    return np.array(rmat_list)


def extract_rlv():
    find_rlv = re.compile(r"recip_lattice\d\s+{(.+)}")
    file_name = open("./Index.txt").read()
    file_split = re.split("[$]pattern\d", file_name)[1:]

    rlv_list = []
    for patt in file_split:
        rlv_string = find_rlv.findall(patt)[0]
        rlv_colum_str = re.findall(r"{[-+]?\d+\.?\d*,[-+]?\d+\.?\d*,[-+]?\d+\.?\d*}", rlv_string)

        rlv_colum_list = []
        for rlv_colum in rlv_colum_str:
            rlv_colum_list.append([float(s) for s in re.findall(r"[-+]?\d+\.?\d*", rlv_colum)])
        rlv_colum_array = np.array(rlv_colum_list)
        rlv_list.append(rlv_colum_array)

    return np.array(rlv_list[0])


# write parameters into a Peaks.txt file
def write_mat_params(intensities, gvectors, material, lattice_params, spacegroup):
    peak_file = open('Peaks.txt', 'w')
    lines = open("trial_gold_Peak.txt").readlines()
    materials = material
    a, b, c, ang1, ang2, ang3 = lattice_params
    space_group = spacegroup
    i = 1
    latticeParameters = "{ %s, %s, %s, %s, %s, %s }" % (a * 1e9, b * 1e9, c * 1e9, ang1, ang2, ang3)
    while i <= 12:
        for s in lines:
            if s.startswith("$filetype") \
                    or s.startswith("// parameters defining the crystal structure") \
                    or s.startswith("$latticeAlphaT") \
                    or s.startswith("$lengthUnit") \
                    or s.startswith("// the following table contains xyz compotnents of G^ and the integral of"
                                    "the peak"):
                peak_file.write(s)
            if s.startswith("$structureDesc"):
                peak_file.write("$structureDesc		%s\n" % materials)
            if s.startswith("$latticeParameters"):
                peak_file.write("$latticeParameters	%s	// 2006, CODATA\n" % latticeParameters)
            if s.startswith("$SpaceGroup"):
                peak_file.write("$SpaceGroup			%s					// Space Group number from International\n"
                                % space_group)
            if s.startswith("$N_Ghat+Intens"):
                peak_file.write("$N_Ghat+Intens 	%s		// number of G^ vectors\n" % len(intensities))
            i += 1
        G_line = 1
        # print("gvectors_group is",gvectors_group)
        for n, new_Value in enumerate(gvectors):
            peak_file.write("%s,%s,%s,%s\n" % (new_Value[0], new_Value[1], new_Value[2], intensities[n]))
        peak_file.close()


def calc_gs(pixel, rot_vec, trans_vec, det_params):
    pix_x, pix_y, pitch_x, pitch_y, name = det_params
    theta = np.sqrt(rot_vec[0] * rot_vec[0] + rot_vec[1] * rot_vec[1] + rot_vec[2] * rot_vec[2])
    c = np.cos(theta)
    s = np.sin(theta)
    c1 = 1 - c
    Rx = rot_vec[0]
    Ry = rot_vec[1]
    Rz = rot_vec[2]
    Rx /= theta
    Ry /= theta
    Rz /= theta
    # print(Rx,Ry,Rz)
    rho00 = c + Rx * Rx * c1
    rho01 = Rx * Ry * c1 - Rz * s
    rho02 = Ry * s + Rx * Rz * c1
    rho10 = Rz * s + Rx * Ry * c1
    rho11 = c + Ry * Ry * c1
    rho12 = -Rx * s + Ry * Rz * c1
    rho20 = -Ry * s + Rx * Rz * c1
    rho21 = Rx * s + Ry * Rz * c1
    rho22 = c + Rz * Rz * c1

    q_vectors = []

    print(rho00, rho01, rho02)
    print(rho10, rho11, rho12)
    print(rho20, rho21, rho22)
    for peak in pixel:
        # peak_xy = self.peak_dict[peak]['XY']
        # start calculating the g-vectors below
        px = peak[0]
        py = peak[1]

        xd = (px - 0.5 * (pix_x - 1)) * pitch_x
        yd = (py - 0.5 * (pix_y - 1)) * pitch_y
        zd = 0
        # translate (xd,yd) by the vector P to get (xd,yd,zd)
        x, y, z = xd + trans_vec[0], yd + trans_vec[1], zd + trans_vec[2]
        # rotate about R
        X, Y, Z = rho00 * x + rho01 * y + rho02 * z, rho10 * x + rho11 * y + rho12 * z, rho20 * x + rho21 * y + rho22 * z
        # normalize
        total = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
        x_lab, y_lab, z_lab = X / total, Y / total, Z / total

        # normalize qhat
        q_x, q_y, q_z = x_lab, y_lab, z_lab - 1
        total_qhat = np.sqrt(q_x ** 2 + q_y ** 2 + q_z ** 2)
        q_x = q_x / total_qhat
        q_y = q_y / total_qhat
        q_z = q_z / total_qhat

        q_vectors.append((q_x, q_y, q_z))

    return q_vectors


def extract_hkl(filename):
    file_name = open(filename).read()

    find_hkls = re.findall(r"\[\s+\d+]\s+\(.+?\)\s+\((.+?)\)", file_name)
    hkls = []
    for hkl in find_hkls:
        hkl = re.findall(r"[\s?](\d+)", hkl)
        hkl_list = []
        for element in hkl:
            element = int(element)
            hkl_list.append(element)
        hkls.append(hkl_list)

    return hkls


# extract_hkl("Index.txt")


def lab_cryst_calc(energy, hkl, gonio_angles):
    h, k, ell = hkl[0], hkl[1], hkl[2]
    delta, theta, chi, phi, gamma = gonio_angles[0], gonio_angles[1], gonio_angles[2], gonio_angles[3], gonio_angles[4]
    # Find kf from beamline code
    kf = ed_qvec(energy, delta, gamma)
    # Normalize
    kf = kf / np.linalg.norm(kf)
    # ki is [0,0,1]
    ki = np.array([0, 0, 1])
    # Calculate Q in lab frame
    q_lab = kf - ki
    # Normalize
    q_lab = q_lab / np.linalg.norm(q_lab)
    # Set Q in crystal frame for the measured reflection
    q_cryst = np.array([h, k, ell]) / np.linalg.norm(np.array([h, k, ell]))

    print('Q_lab:', q_lab, '\nQ_cryst:', q_cryst)

    # # Rotate Q_lab into Crystal Frame via axis-angle rotation
    # ax = np.cross(q_lab, q_cryst)  # find axis by taking the cross product of Q_lab and Q_cryst
    # ax = ax / np.linalg.norm(ax)
    # ang = np.arccos(q_cryst @ q_lab)  # find angle by taking arccos of dot product of Q_lab and Q_cryst
    # o_mat = Rotation.from_rotvec(ax * ang).as_matrix()  # convert to orientation matrix
    # # create list of ax-angle rotations about Q_cryst
    # ax_ang = [ang * q_cryst for ang in np.linspace(np.pi / 4, np.pi / 2, 100)]
    # # convert ax-angles to rotation matrices and dot with our orientaion mat
    # o_list = [Rotation.from_rotvec(r).as_matrix() @ o_mat for r in ax_ang]
    # o_list = [o.T for o in o_list]  # Transpose all orientations

    return np.array(q_lab), np.array(q_cryst)


def ed_qvec(energy, delta, gamma):
    deg2rad = np.pi / 180
    lam = (12.398 / energy) / 10
    delta = delta * deg2rad
    gam = gamma * deg2rad
    shift1, shift2, shift3 = 500.0, 500.0, 500.0

    Qlabcenter1 = np.sin(delta) * np.cos(gam) * (2 * 3.14159265) / lam
    Qlabcenter2 = np.sin(gam) * (2 * 3.14159265) / lam
    Qlabcenter3 = (np.cos(delta) * np.cos(gam) - 1.0) * (2 * 3.14159265) / lam
    ki = (0.0, 0.0, 1.0)
    kf = (np.sin(delta) * np.cos(gam) * (2 * 3.14159265) / lam, np.sin(gam) * (2 * 3.14159265) / lam,
          np.cos(delta) * np.cos(gam) * (2 * 3.14159265) / lam)
    # print (ki,kf)
    # print (Qlabcenter1, Qlabcenter2, Qlabcenter3)
    # vectorarr = vtk.vtkDoubleArray()
    # vectorarr.SetNumberOfComponents(3)
    # vectorarr.SetNumberOfTuples(3)
    # vectorarr.SetComponent(0,0,Qlabcenter1)
    # vectorarr.SetComponent(0,1,Qlabcenter2)
    # vectorarr.SetComponent(0,2,Qlabcenter3)
    # vectorarr.SetComponent(1,0,ki[0])
    # vectorarr.SetComponent(1,1,ki[1])
    # vectorarr.SetComponent(1,2,ki[2])
    # vectorarr.SetComponent(2,0,kf[0])
    # vectorarr.SetComponent(2,1,kf[1])
    # vectorarr.SetComponent(2,2,kf[2])
    #
    # print (vectorarr)
    # vectorpoints=vtk.vtkPoints()
    # vectorpoints.SetDataTypeToDouble()
    # vectorpoints.SetNumberOfPoints(3)
    # vectorpoints.SetPoint(0, (shift1,shift2,shift3) )
    # vectorpoints.SetPoint(1, (shift1,shift2,shift3) )
    # vectorpoints.SetPoint(2, (shift1,shift2,shift3) )
    # print (vectorpoints)
    #
    # vectorgrid = vtk.vtkUnstructuredGrid()
    # vectorgrid.SetPoints(vectorpoints)
    # vectorgrid.GetPointData().SetVectors(vectorarr)
    #
    # usgridwriter=vtk.vtkUnstructuredGridWriter()
    # usgridwriter.SetFileName('Qvector.vtk')
    # usgridwriter.SetFileTypeToASCII()
    # usgridwriter.SetInputData(vectorgrid)
    # usgridwriter.Write()
    return kf

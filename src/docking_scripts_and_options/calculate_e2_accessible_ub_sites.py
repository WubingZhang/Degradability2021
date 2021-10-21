__author__ = "Shourya S. Roy Burman"
__email__ = "ssrburman@gmail.com"

import os
import argparse
import numpy as np
import pandas as pd
import time
from pymol import cmd, util


def align_ref(model, rad):
    """
    Return distance between CRBN-LVY N3 and aligned reference 2PU C26
    :param model: Name of pymol selection
    :param rad: Radius within which free path has to be calculated
    :return: True if alignment works, False otherwise
    """
    cmd.super("ref and resi 1-149", "{} and chain B".format(model))
    dist = cmd.distance("{} and chain A and resn LVY and name N4".format(model), "ref and resn GUI and name N01")
    num_bb = has_free_path(model, rad)
    return dist, num_bb


def is_in_cylinder(a, b, rad, p):
    """
    Returns True if point p is in cylinder of radius rad joining points a and b
    Math from https://bit.ly/3oRNAdN
    :param a: numpy array coordinates of a
    :param b: numpy array coordinates of b
    :param rad: radius of cylinder
    :param p: numpy array coordinates of p
    :return: Boolean
    """

    e = b - a
    m = np.cross(a, b)
    d = np.linalg.norm(m + np.cross(e, p)) / np.linalg.norm(e)  # distance of p from line joining a and b
    if d >= rad:
        return False
    else:
        q = p + (np.cross(e, m + np.cross(e, p))) / (np.linalg.norm(e)) ** 2  # q is closest point to p on line
        # calculating barycentric coordinates of q
        wa = np.linalg.norm(np.cross(q, b)) / np.linalg.norm(m)
        wb = np.linalg.norm(np.cross(q, a)) / np.linalg.norm(m)
        # is q between a and b?
        if 0 <= wa <= 1 and 0 <= wb <= 1:
            return True
        else:
            return False


def has_free_path(model, rad):
    """
    Determine if there is a free path between the two points. Only protein backbone atoms and cmpd non-neighbours are
    considered
    :param model: Name of pymol selection
    :param rad: Radius if the cylinder within which to calculate clashes
    :return: Number of clashing atoms in cylinder
    """

    a = cmd.get_coords("{} and chain A and resn LVY and name N4".format(model), 1)[0]
    b = cmd.get_coords("ref and resn GUI and name N01", 1)[0]
    bb_coords = cmd.get_coords("{} and name CA+C+N".format(model), 1)
    lvy_coords = cmd.get_coords("{} and resn LVY and !(name N4+C4+C3+C5)".format(model))  # LVY except neighboring atoms
    gui_coords = cmd.get_coords("ref and resn GUI and !(name N01+C18+C20)")  # GUI except neighboring atoms
    num_bb = 0
    for bb in np.concatenate([bb_coords, lvy_coords, gui_coords]):
        if is_in_cylinder(a, b, rad, bb):
            num_bb += 1

    return num_bb


def get_plane_coeffs(model):
    """
    Calculate vertical plane passing through 1/2 dist between COM_target and CRBN, and hortizontal plane passing through
    COM_target
    :param model: Name of pymol selection
    :return: A pair of 4-tupples for the coefficients for each plane
    """

    # Vertical: Getting coordinates for C-alpha of CRBN R130, Y151, and I333
    ver_a = np.array(cmd.get_coords(selection="{} and chain A and resi 130 and name CA".format(model))[0])
    ver_b = np.array(cmd.get_coords(selection="{} and chain A and resi 151 and name CA".format(model))[0])
    ver_c = np.array(cmd.get_coords(selection="{} and chain A and resi 333 and name CA".format(model))[0])

    # Pass planes through centre of mass of kinase
    com_kinase = np.array(cmd.centerofmass("{} and chain B".format(model)))

    # These two vectors are in the plane
    v1 = ver_c - ver_a
    v2 = ver_b - ver_a

    # the cross product is a vector normal to the plane
    ver_cp = np.cross(v1, v2)
    a, b, c = ver_cp

    # This evaluates a * x_c + b * y_c + c * z_c which equals d
    # d_on_plane = np.dot(ver_cp, ver_c)
    # Taking projection of COM on vertical plane passing through CRBN
    t = (np.dot(ver_cp, ver_c) - np.dot(ver_cp, com_kinase)) / np.linalg.norm(ver_cp) ** 2
    com_proj = np.array([com_kinase[0] + t * a, com_kinase[1] + t * b, com_kinase[2] + t * c])
    com_half_dist = np.mean(np.array([com_kinase, com_proj]), axis=0)
    d = np.dot(ver_cp, com_half_dist)
    # Taking plane through half distance between CRBN and COM
    ver_plane_coeffs = np.array([a, b, c, -d])

    # print("The vertical plane is {}x + {}y + {}z + {} = 0".format(*ver_plane_coeffs))

    # Horizontal: Getting coordinates for C-alpha of CRBN D64, 208, 409
    hor_a = np.array(cmd.get_coords(selection="{} and chain A and resi 64 and name CA".format(model))[0])
    hor_b = np.array(cmd.get_coords(selection="{} and chain A and resi 208 and name CA".format(model))[0])
    hor_c = np.array(cmd.get_coords(selection="{} and chain A and resi 409 and name CA".format(model))[0])

    # These two vectors are in the plane
    v3 = hor_c - hor_a
    v4 = hor_b - hor_a

    # the cross product is a vector normal to the plane
    hor_cp = np.cross(v3, v4)
    p, q, r = hor_cp

    # This evaluates p * x_c + q * y_c + r * z_c which equals s
    s = np.dot(hor_cp, com_kinase)
    hor_plane_coeffs = np.array([p, q, r, -s])

    # print("The horizontal plane is {}x + {}y + {}z + {} = 0".format(*hor_plane_coeffs))

    return ver_plane_coeffs, hor_plane_coeffs


def is_e2_accessible(p1, ver_plane_coeffs, hor_plane_coeffs, lvy_n4):
    """
    Is a point e2 acccessible?
    :param p1: 3-D coordinats of point for which the calculations are being made
    :param ver_plane_coeffs: Coefficients of vertical plane
    :param hor_plane_coeffs:Coefficients of horizontal plane
    :param lvy_n4: 3-D coordinates of N4 atom of LVY
    :return: Boolean
    """

    if round(np.dot(ver_plane_coeffs, np.append(p1, 1))) <= 0 <= round(np.dot(hor_plane_coeffs, np.append(p1, 1))) and \
        np.linalg.norm(lvy_n4 - p1) <= 60:  # Rounding to take care of numerical errors
        return True
    else:
        return False


def lys_ca_cb_orientation(p0, p1, ver_plane_coeffs, hor_plane_coeffs):
    """
    Is the lysine Ca->Cb vector pointing in the completely wrong direction?
    :param p0: 3-D coordinates of C-alpha
    :param p1: 3-D coordinates of C-beta
    :param ver_plane_coeffs: Coefficients of vertical plane
    :param hor_plane_coeffs: Coefficients of horizontal plane
    :return:
    """

    ver_orient = np.dot(ver_plane_coeffs[:3], p1 - p0)
    hor_orient = np.dot(hor_plane_coeffs[:3], p1 - p0)

    if hor_orient < 0 < ver_orient:  # Both orientations are wrong
        return False
    else:
        return True


def get_ub_sites(ub_sites_df, offset_df, src_dir, dst):
    """
    Writes file with fraction of E2 accessibility of every Ub site in the ub_sites dataframe.
    Also, writes .pse files highlighting which Ub sites are considered accessible, what which ones are not.
    :param ub_sites_df: Input dataframe of Ub sites
    :param offset_df: PDB numbering offsets (pre-calculated)
    :param src_dir: Directory where all files are kept
    :param dst: Directory where
    """
    ub_sites = {}
    pdb_with_offsets = set()

    with open(os.path.join(src_dir, "e2_accessible_sites.csv"), "w") as f:
        f.write("Uniprot_AC,Ub_site,E2_accessibility,E2_accessible_orientation\n")

        for pdb in sorted(os.listdir(src_dir)):
            if not pdb.startswith("."):
                start = time.time()
                ac = pdb.split("_")[0]
                if offset_df[offset_df["PDB"] == pdb]["Offset"].count() == 0:  # No offset
                    offset = 0
                else:
                    offset = offset_df[offset_df["PDB"] == pdb]["Offset"].to_list()[0]
                ub_sites[ac] = [int(site.split("-")[0][1:]) + offset
                                for site in ub_sites_df[ub_sites_df["ACC_ID"] == ac]["MOD_RSD"].to_list()
                                if int(site.split("-")[0][1:]) + offset > 0]
                if not os.path.exists(os.path.join(src_dir, pdb, dst)):
                    os.mkdir(os.path.join(src_dir, pdb, dst))

                # Count accessible site
                e2_accessible_site_count = dict.fromkeys(ub_sites[ac], None)  # blank dictionary with Ub sites
                e2_accessible_site_orient_count = dict.fromkeys(ub_sites[ac], None)  # blank dictionary with Ub sites
                model_count = 0  # PROTAC feasible model count

                for model in os.listdir(os.path.join(src_dir, pdb, "top_models")):
                    cmd.load(os.path.join(src_dir, pdb, "top_models", model), model[:-7])
                    util.cbc(model[:-7])

                    # Is there a free path from LVY? Max dist: 14 Ang, max number of atoms in cylinder:2, radius: 1 Ang
                    dist, num_bb = align_ref(model[:-7], 1)

                    if dist <= 14 and num_bb <= 2:  # PROTAC feasible

                        model_count += 1
                        ver_plane_coeffs, hor_plane_coeffs = get_plane_coeffs(model[:-7])
                        lvy_n4 = np.array(cmd.get_coords(
                            selection="{} and chain A and resn LVY and name N4".format(model[:-7]))[0])

                        # Check that they are indeed lysine
                        for site in ub_sites[ac]:
                            resi_names = {"resi_name": []}
                            cmd.iterate("{} and chain B and resi {} and name CA".format(model[:-7], site),
                                        "resi_name.append(resn)", space=resi_names)
                            if len(resi_names["resi_name"]) > 0:  # residue might not exist in structure
                                if resi_names["resi_name"][0] != "LYS":
                                    print("In {}, at {}, residue in {}.".format(pdb, site, resi_names["resi_name"][0]))
                                    pdb_with_offsets.add(pdb)
                                else:  # If they are lysine
                                    p0 = np.array(
                                        cmd.get_coords(selection="{} and chain B and resi {} and name CA".format(
                                            model[:-7], site))[0])
                                    p1 = np.array(
                                        cmd.get_coords(selection="{} and chain B and resi {} and name CB".format(
                                            model[:-7], site))[0])

                                    # Color ub sites by E2 accessibility
                                    if is_e2_accessible(p1, ver_plane_coeffs, hor_plane_coeffs, lvy_n4):
                                        cmd.show("sticks", "{} and chain B and resi {}".format(model[:-7], site))
                                        cmd.color("blue", "{} and chain B and resi {}".format(model[:-7], site))
                                        if e2_accessible_site_count[site] is None:
                                            e2_accessible_site_count[site] = 1
                                        else:
                                            e2_accessible_site_count[site] += 1

                                        if lys_ca_cb_orientation(p0, p1, ver_plane_coeffs, hor_plane_coeffs):
                                            if e2_accessible_site_orient_count[site] is None:
                                                e2_accessible_site_orient_count[site] = 1
                                            else:
                                                e2_accessible_site_orient_count[site] += 1

                                    else:
                                        cmd.show("sticks", "{} and chain B and resi {}".format(model[:-7], site))
                                        cmd.color("firebrick", "{} and chain B and resi {}".format(model[:-7], site))
                                        if e2_accessible_site_count[site] is None:
                                            e2_accessible_site_count[site] = 0
                                            e2_accessible_site_orient_count[site] = 0

                        util.cnc(model[:-7])
                        cmd.save(os.path.join(src_dir, pdb, dst, "{}_e2_accessible_sites.pse".format(model[:-7])),
                                 "{} + ref".format(model[:-7]))
                    cmd.delete(model[:-7])

                end = time.time()
                print("It took {} seconds to finish {}.".format(end - start, pdb))

                for site, count in e2_accessible_site_count.items():
                    if count is not None and model_count > 0:
                        if e2_accessible_site_orient_count[site] is not None:
                            f.write("{},{},{},{}\n".format(pdb.split("_")[0], site, (count / model_count),
                                                           (e2_accessible_site_orient_count[site] / model_count)))
                        else:
                            f.write("{},{},{},{}\n".format(pdb.split("_")[0], site, (count / model_count),
                                                           e2_accessible_site_orient_count[site]))
                    else:
                        f.write("{},{},{},{}\n".format(pdb.split("_")[0], site, count,
                                                       e2_accessible_site_orient_count[site]))
                f.flush()


def main():
    # Parsing all arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--src_dir", help="Super directory with clustered decoys",
                        default="../07_local_docking_all_kinases/01_CRBN_Len")
    parser.add_argument("-d", "--dst", help="Where the pymol sessions will be stored",
                        default="e2_accessible_ub_sites")
    parser.add_argument("-u", "--ub_sites", help="File with ubiquitination sites",
                        default="Ubiquitination_site_dataset")
    parser.add_argument("-r", "--ref_struct", help="Name of reference kinase to align",
                        default="../07_local_docking_all_kinases/aux_files/CDK2_TAE.pdb")
    parser.add_argument("-o", "--offset_file", help="Helper file with offsets for special cases",
                        default="PDB_numbering_offsets.csv")

    args = parser.parse_args()
    src_dir = args.src_dir
    ub_sites_df = pd.read_csv(args.ub_sites, delimiter="\t", header=2, low_memory=False)
    offset_df = pd.read_csv(args.offset_file)
    dst = args.dst

    # Loading fixed reference
    cmd.load(os.path.join(args.ref_struct), "ref")
    cmd.color("gray90", "ref")
    cmd.set("cartoon_transparency", 0.4, "ref")
    util.cnc("ref")

    get_ub_sites(ub_sites_df, offset_df, src_dir, dst)


main()

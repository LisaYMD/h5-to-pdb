import sys
import os
import numpy as np
from tqdm import tqdm
import json
import h5py

class h5toPDB:
    def __init__(self, filename):
        # load h5 file
        f = h5py.File(filename, 'r')
        self.filename = filename
        # read types and constants
        types = f['readdy']['config']['particle_types']
        self.moltype = [types[k][0].decode() for k in range(0, len(types))] 
        self.molnum = [types[k][1] for k in range(0, len(types))]
        self.mol_diffusion = [types[k][2] for k in range(0, len(types))]
        # read trajectory
        if 'trajectory' not in f['readdy']:
            raise ValueError("No trajectory found in the provided h5 file.")
        self.limrec = f['readdy']['trajectory']['limits']
        self.rec = f['readdy']['trajectory']['records']
        self.traj_len = int(self.limrec.shape[0]-1)
        # read topologies
        if 'topologies' not in f['readdy']['observables']:
            raise ValueError("No topologies found in the provided h5 file.")
        topo = f['readdy']['observables']['topologies']
        self.limparts = topo['limitsParticles']
        self.limedges = topo['limitsEdges']
        self.parts = topo['particles'] 
        self.edges = topo['edges']
        # read box size
        info = json.loads(np.atleast_1d(f['readdy']['config']['general'][...])[0].decode()) # box_size, box_volume, kbt, pbc
        self.kbt = info["kbt"]
        self.xbox, self.ybox, self.zbox = info["box_size"][0], info["box_size"][1], info["box_size"][2]
        # molecule size and color
        eta = 0.85137 # assume the viscosity of water (milliPascal*second) # "cell"
        self.molsize = [self.kbt*10/(6.02214076*6*np.pi*eta*self.mol_diffusion[k]) for k in range(0, len(self.mol_diffusion))]
        self.molcolor = [0 for k in range(0, len(self.moltype))] # default color = 0 (blue) 
        return None

    def connected_edge(self, tim):
        timedges = self.edges[self.limedges[tim,0]:self.limedges[tim,1]]
        timparts = self.parts[self.limparts[tim,0]:self.limparts[tim,1]]
        index = 0
        # first, change local edgelists into global edgelists
        local_partlist = []
        global_edgelist = np.empty((0,2))
        while index < (self.limparts[tim,1]-self.limparts[tim,0]):
            topvalue = timparts[index]
            local_partlist.append(np.array(timparts[int(index+1):int(index+topvalue+1)]))
            index = int(index+topvalue+1)
        index2 = 0
        i = 0
        while index2 < (self.limedges[tim,1]-self.limedges[tim,0]):
            edgevalue = timedges[index2,0]
            local_edge = timedges[int(index2+1):int(index2+edgevalue+1)]
            global_part = local_partlist[i][:]
            lfunc = lambda x: global_part[x]
            global_part = local_partlist[i][:]
            global_edge = lfunc(local_edge)
            global_edgelist = np.concatenate([global_edgelist, global_edge])
            index2 = int(index2+edgevalue+1)
            i += 1
        return global_edgelist
 
    ## detect only existing types and molecules
    # update moltype, molnum, molsize, molcolor
    def update_existingtypes(self):
        accum_types = set()
        for t in range(0, self.traj_len):
            current_types = set([self.rec[m][0] for m in range(self.limrec[t,0], self.limrec[t,1])])
            accum_types = accum_types | current_types
        molnum2 = list(accum_types)
        moltype2 = [self.moltype[self.molnum.index(k)] for k in molnum2]
        molsize2 = [self.molsize[self.molnum.index(k)] for k in molnum2]
        molcolor2 = [self.molcolor[self.molnum.index(k)] for k in molnum2]
        return moltype2, molnum2, molsize2, molcolor2

    ## if you want to paint different molecule as different colors, you can consult setting files
    def resid_from_setting(self):
        if os.path.exists("mol_settings.json"):
            with open("mol_settings.json") as json_file:
                settings = json.load(json_file)
            resid_candidates = settings["molecule"]
            molcompose = settings["particle_type"]
            # sort by moletypes
            moltype_settings, molsize_update, molcolor_update = [], [], []
            for r in resid_candidates:
                moltype_settings.extend(settings["particle_type"][r])
                molsize_update.extend(settings["size"][r])
                molcolor_update.extend(settings["color"][r])
            molsize = [molsize_update[self.moltype.index(m)] for m in moltype_settings if m in self.moltype]
            molcolor = [molcolor_update[self.moltype.index(m)] for m in moltype_settings if m in self.moltype]
        else:
            resid_candidates = []
            molcompose = []
            molsize = self.molsize
            molcolor = self.molcolor
        return resid_candidates, molcompose, molsize, molcolor
    
    def generate_pdb(self):
        print("Generating PDB file...")
        start, duration = 0, self.traj_len
        edge_data = self.connected_edge(0)
        self.moltype, self.molnum, self.molsize, self.molcolor = self.update_existingtypes()
        resid_candidates, molcompose, _, _ = self.resid_from_setting()
        fname = self.filename[:-3] + ".pdb"
        # record the writing strings in list "lines"
        resid = "PSD" # by default
        ### Generate PDB file
        with open(fname, "w") as f:
            f.write("CRYST1  "+f"{format(self.xbox*10, '.2f'):<8}"+" "+f"{format(self.ybox*10, '.2f'):<8}"+" "+f"{format(self.zbox*10, '.2f'):<8}"+" 90.00  90.00  90.00\n") 
            for t in tqdm(range(start, start+duration)):
                lines = []
                lines.append("MODEL\n")
                for m in range(self.limrec[t,0], self.limrec[t,1]):
                    pos0 = str(format(10*self.rec[m][3][0], '.1f')) 
                    pos1 = str(format(10*self.rec[m][3][1], '.1f')) 
                    pos2 = str(format(10*self.rec[m][3][2], '.1f')) 
                    position = " "*(7-len(pos0)+4)+ pos0 + " "*(7-len(pos1)+1) + pos1 + " "*(7-len(pos2)+1) + pos2
                    if self.rec[m][1] < 10:
                        space = "   "
                    elif self.rec[m][1] < 100:
                        space = "  "
                    elif self.rec[m][1] < 1000:
                        space = " "
                    else:
                        space = ""
                    if resid_candidates != []:
                        for r in range(0, len(resid_candidates)):
                            if self.moltype[self.molnum.index(self.rec[m][0])] in molcompose[resid_candidates[r]]:
                                resid = resid_candidates[r][:3]
                                break
                    lines.append("ATOM   "+space+str(self.rec[m][1])+" "+f"{self.moltype[self.molnum.index(self.rec[m][0])][:4]:<4}"+" "+resid+" A"+space+str(self.rec[m][1])+position+"  0.00  0.00\n")
                lines.append("ENDMDL\n")
                f.write("".join(lines))
            f.write("END")
        print("PDB file generated: " + fname)
        return None

    ### necessary for drawing edge
    def generate_tcl(self):
        print("Generating Tcl script file...")
        resid_candidates, molcompose, _, _ = self.resid_from_setting()
        self.moltype, self.molnum, self.molsize, self.molcolor = self.update_existingtypes() 
        edge_list = []
        edge_data = self.connected_edge(0)
        for e in range(0, len(edge_data)):
            edge_list.append(tuple(edge_data[e,:].astype(int))) 
        ### Generate Tcl script file
        fname_tcl = self.filename[:-3] + ".pdb.tcl"
        with open(fname_tcl, 'w') as g:
            lines_tcl = []
            lines_tcl.append("mol delete top\n")
            lines_tcl.append("mol load pdb "+fname_tcl[:-4]+"\n")
            lines_tcl.append("mol delrep 0 top\n")
            lines_tcl.append("display resetview\n")
            overlap = np.zeros(len(self.moltype))
            for l in range(0, len(self.moltype)):
                lines_tcl.append("mol representation VDW "+str(self.molsize[l]*6.6)+"\n")
                lines_tcl.append("mol selection name "+self.moltype[l]+"\n")
                if resid_candidates != []:
                    for r in range(0, len(resid_candidates)):
                        if self.moltype[l] in molcompose[resid_candidates[r]]:
                            lines_tcl.append("mol color ColorID "+str(int(self.molcolor[l]))+"\n")
                            break
                else:
                    lines_tcl.append("mol color ColorID "+str(int(self.molcolor[l]))+"\n")
                lines_tcl.append("mol addrep top\n")
            lines_tcl.append("animate goto 0\n")
            lines_tcl.append("color Display Background white\n")
            lines_tcl.append("molinfo top set {center_matrix} {{{1 0 0 0}{0 1 0 0}{0 0 1 0}{0 0 0 1}}}\n")
            lines_tcl.append("set x "+str(self.xbox*10)+"\n")
            lines_tcl.append("set y "+str(self.ybox*10)+"\n")
            lines_tcl.append("set z "+str(self.zbox*10)+"\n")
            # draw lines between each connected domains
            # adding bonds
            lines_tcl.append("set sel [atomselect $top all]\n")
            for e in range(0, len(edge_list)):
                lines_tcl.append("topo addbond "+str(edge_list[e][1])+" "+str(edge_list[e][0])+"\n")
            # drawing bonds
            lines_tcl.append("mol representation Bonds 1.5\n")
            if resid_candidates != []:
                for r in range(0, len(resid_candidates)):
                    lines_tcl.append("mol selection resname "+str(resid_candidates[r][:3])+"\n")
                    lines_tcl.append("mol color ColorID "+str(int(r+1))+"\n")
                    lines_tcl.append("mol addrep top\n")
            else:
                lines_tcl.append("mol selection resname PSD\n")
                lines_tcl.append("mol color ColorID 2\n")
                lines_tcl.append("mol addrep top\n")
            # drawing pbc box
            lines_tcl.append("pbc set {$x $y $z 90.0 90.0 90.0 -all}\n")
            lines_tcl.append("pbc box_draw -center origin -color black\n")
            # if necesssary: wrap bonds
            g.write("".join(lines_tcl))
        print("Tcl script file generated: " + fname_tcl)
        return None
    
    ## generate lammps data file
    def generate_lammpstrj(self):
        print("Generating LAMMPS trajectory file...")
        start, duration = 0, self.traj_len
        edge_data = self.connected_edge(0)
        self.moltype, self.molnum, self.molsize, self.molcolor = self.update_existingtypes()
        resid_candidates, molcompose, self.molsize, self.molcolor = self.resid_from_setting()
        fname = self.filename[:-3] + ".lammpstrj"
        with open(fname, "w") as f:
            for t in tqdm(range(start, start+duration)):
                lines = []
                lines.append("ITEM: TIMESTEP\n"+str(t)+"\nITEM: NUMBER OF ATOMS\n"+str(self.limrec[t,1]-self.limrec[t,0])+"\nITEM: BOX BOUNDS pp pp pp\n")
                lines.append(str(-self.xbox/2)+" "+str(self.xbox/2)+"\n"+str(-self.ybox/2)+" "+str(self.ybox/2)+"\n"+str(-self.zbox/2)+" "+str(self.zbox/2)+"\n")
                lines.append("ITEM: ATOMS id type xs ys zs\n")
                for m in range(self.limrec[t,0], self.limrec[t,1]):
                    pos0 = str(format(self.rec[m][3][0], '.6f'))
                    pos1 = str(format(self.rec[m][3][1], '.6f'))
                    pos2 = str(format(self.rec[m][3][2], '.6f'))
                    resid = 1
                    lines.append(str(self.rec[m][1])+" "+str(resid)+" "+pos0+" "+pos1+" "+pos2+"\n")
                f.write("".join(lines))
        print("LAMMPS trajectory file generated: " + fname)
        return None

if __name__ == '__main__':
    filename = sys.argv[1]
    result = h5toPDB(filename)
    result.generate_pdb()
    result.generate_tcl()

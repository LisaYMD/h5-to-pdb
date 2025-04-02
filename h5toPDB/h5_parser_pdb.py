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
        self.limrec = f['readdy']['trajectory']['limits']
        self.rec = f['readdy']['trajectory']['records']
        self.traj_len = int(self.limrec.shape[0]-1)
        # read topologies
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
            molsize_update = self.molsize
            molcolor_update = self.molcolor
        return resid_candidates, molcompose, molsize, molcolor
    
    def generate_pdb(self):
        start, duration = 0, self.traj_len
        edge_data = self.connected_edge(0)
        self.moltype, self.molnum, self.molsize, self.molcolor = self.update_existingtypes()
        resid_candidates, molcompose, _, _ = self.resid_from_setting()
        fname = self.filename[:-3] + ".pdb"
        f = open(fname, "w")
        f.write("CRYST1"+"  "+f"{format(self.xbox*10, '.3f'):<8}"+" "+f"{format(self.ybox*10, '.3f'):<8}"+" "+f"{format(self.zbox*10, '.3f'):<8}"+" 90.00  90.00  90.00"+"\n")
        ### Generate PDB file
        for t in tqdm(range(start, start+duration)):
            f.write("MODEL\n")
            for m in range(self.limrec[t,0], self.limrec[t,1]):
                pos0 = str(format(10*self.rec[m][3][0], '.1f')) # 166.2
                pos1 = str(format(10*self.rec[m][3][1], '.1f')) # 218.1
                pos2 = str(format(10*self.rec[m][3][2], '.1f')) # -412.2
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
                else:
                    resid = "PSD"
                f.write("ATOM   "+space+str(self.rec[m][1])+" "+f"{self.moltype[self.molnum.index(self.rec[m][0])][:4]:<4}"+" "+resid+" A"+space+str(self.rec[m][1])+position+"  0.00  0.00\n")
            f.write("ENDMDL\n")
        f.write("END") 
        f.close()
        return None

    ### necessary for drawing edge
    def generate_tcl(self):
        resid_candidates, molcompose, self.molsize, self.molcolor = self.resid_from_setting() 
        edge_list = []
        edge_data = self.connected_edge(0)
        for e in range(0, len(edge_data)):
            edge_list.append(tuple(edge_data[e,:].astype(int))) 
        ### Generate Tcl script file
        fname = self.filename[:-3] + ".pdb"
        fname_tcl = self.filename[:-3] + ".pdb.tcl"
        with open(fname_tcl, 'w') as g:
            g.write("mol delete top\n")
            g.write("mol load pdb "+fname+"\n")
            g.write("mol delrep 0 top\n")
            g.write("display resetview\n")
            overlap = np.zeros(len(self.moltype))
            for l in range(0, len(self.moltype)):
                g.write("mol representation VDW "+str(self.molsize[l]*6.6)+"\n")
                g.write("mol selection name "+self.moltype[l]+"\n")
                if resid_candidates != []:
                    for r in range(0, len(resid_candidates)):
                        if self.moltype[l] in molcompose[resid_candidates[r]]:
                            g.write("mol color ColorID "+str(int(self.molcolor[l]))+"\n")
                            break
                else:
                    g.write("mol color ColorID "+str(int(self.molcolor[l]))+"\n")
                g.write("mol addrep top\n")
            g.write("animate goto 0\n")
            g.write("color Display Background white\n")
            g.write("molinfo top set {center_matrix} {{{1 0 0 0}{0 1 0 0}{0 0 1 0}{0 0 0 1}}}\n")
            g.write("set x "+str(self.xbox*10)+"\n")
            g.write("set y "+str(self.ybox*10)+"\n")
            g.write("set z "+str(self.zbox*10)+"\n")
            # draw lines between each connected domains
            # adding bonds
            g.write("set sel [atomselect $top all]\n")
            for e in range(0, len(edge_list)):
                g.write("topo addbond "+str(edge_list[e][1])+" "+str(edge_list[e][0])+"\n")
            # drawing bonds
            g.write("mol representation Bonds 1.5\n")
            if resid_candidates != []:
                for r in range(0, len(resid_candidates)):
                    g.write("mol selection resname "+str(resid_candidates[r][:3])+"\n")
                    g.write("mol color ColorID "+str(int(r+1))+"\n")
                    g.write("mol addrep top\n")
            else:
                g.write("mol selection resname PSD\n")
                g.write("mol color ColorID 2\n")
                g.write("mol addrep top\n")
            # drawing pbc box
            g.write("pbc set {$x $y $z 90.0 90.0 90.0 -all}\n")
            g.write("pbc box_draw -center origin -color black\n")
            # if necesssary: wrap bonds
        return None
    
    ## generate lammps data file
    def generate_lammpstrj(self):
        start, duration = 0, self.traj_len
        edge_data = self.connected_edge(0)
        self.moltype, self.molnum, self.molsize, self.molcolor = self.update_existingtypes()
        resid_candidates, molcompose, self.molsize, self.molcolor = self.resid_from_setting()
        fname = self.filename[:-3] + ".lammpstrj"
        f = open(fname, "w")
        for t in tqdm(range(start, start+duration)):
            f.write("ITEM: TIMESTEP\n"+str(t)+"\nITEM: NUMBER OF ATOMS\n"+str(self.limrec[t,1]-self.limrec[t,0])+"\nITEM: BOX BOUNDS pp pp pp\n")
            f.write(str(-self.xbox/2)+" "+str(self.xbox/2)+"\n"+str(-self.ybox/2)+" "+str(self.ybox/2)+"\n"+str(-self.zbox/2)+" "+str(self.zbox/2)+"\n")
            f.write("ITEM: ATOMS id type xs ys zs\n")
            for m in range(self.limrec[t,0], self.limrec[t,1]):
                pos0 = str(format(self.rec[m][3][0], '.6f'))
                pos1 = str(format(self.rec[m][3][1], '.6f'))
                pos2 = str(format(self.rec[m][3][2], '.6f'))
                resid = 1
                f.write(str(self.rec[m][1])+" "+str(resid)+" "+pos0+" "+pos1+" "+pos2+"\n")
        f.close()
        return None

if __name__ == '__main__':
    filename = sys.argv[1]
    result = h5toPDB(filename)
    result.generate_pdb()
    result.generate_tcl()

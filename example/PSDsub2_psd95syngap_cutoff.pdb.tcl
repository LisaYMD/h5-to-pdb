mol delete top
mol load pdb PSDsub2_psd95syngap_cutoff.pdb
mol delrep 0 top
display resetview
mol representation VDW 2.7786
mol selection name C
mol color ColorID 10
mol addrep top
mol representation VDW 8.4744
mol selection name D01
mol color ColorID 22
mol addrep top
mol representation VDW 8.4744
mol selection name D02
mol color ColorID 22
mol addrep top
mol representation VDW 8.4744
mol selection name D03
mol color ColorID 22
mol addrep top
mol representation VDW 8.4744
mol selection name D11
mol color ColorID 22
mol addrep top
mol representation VDW 8.4744
mol selection name D12
mol color ColorID 22
mol addrep top
mol representation VDW 8.4744
mol selection name D13
mol color ColorID 22
mol addrep top
mol representation VDW 7.814399999999999
mol selection name E
mol color ColorID 23
mol addrep top
mol representation VDW 10.942799999999998
mol selection name F0
mol color ColorID 0
mol addrep top
mol representation VDW 10.942799999999998
mol selection name PH
mol color ColorID 0
mol addrep top
mol representation VDW 16.6518
mol selection name C2
mol color ColorID 32
mol addrep top
mol representation VDW 12.9954
mol selection name GAP
mol color ColorID 32
mol addrep top
mol representation VDW 13.2
mol selection name Coiltip_sngp
mol color ColorID 32
mol addrep top
mol representation VDW 3.3
mol selection name PBM0
mol color ColorID 4
mol addrep top
mol representation VDW 6.6
mol selection name PBM1
mol color ColorID 32
mol addrep top
animate goto 0
color Display Background white
molinfo top set {center_matrix} {{{1 0 0 0}{0 1 0 0}{0 0 1 0}{0 0 0 1}}}
set x 484.3658212699539
set y 484.3658212699539
set z 484.3658212699539
set sel [atomselect $top all]
topo addbond 1 0
topo addbond 2 1
topo addbond 3 2
topo addbond 4 3
topo addbond 5 4
topo addbond 7 6
topo addbond 8 7
topo addbond 9 8
topo addbond 10 9
topo addbond 11 10
topo addbond 13 12
topo addbond 14 13
topo addbond 15 14
topo addbond 16 15
topo addbond 17 16
topo addbond 28 27
topo addbond 19 18
topo addbond 20 19
topo addbond 22 21
topo addbond 23 22
topo addbond 25 24
topo addbond 26 25
topo addbond 27 20
topo addbond 27 23
topo addbond 27 26
topo addbond 29 28
topo addbond 30 28
topo addbond 31 28
mol representation Bonds 1.5
mol selection resname psd
mol color ColorID 1
mol addrep top
mol selection resname syn
mol color ColorID 2
mol addrep top
pbc set {$x $y $z 90.0 90.0 90.0 -all}
pbc box_draw -center origin -color black

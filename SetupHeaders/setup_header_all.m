% Data names and colors.
nameNMBC={'wt','epe1','swi6','gcn5','set2','set1'};
nameMBC={'wt_MBC','epe1_MBC','swi6_MBC','gcn5_MBC','set2_MBC','set1_MBC'};
nameIntercleave=reshape([nameNMBC;nameMBC],1,2*length(nameNMBC));
nameAll=[nameNMBC,nameMBC];
nameAll2 = {'WT','\itepe1\Delta','\itswi6\Delta','\itgcn5\Delta',...
    '\itset2\Delta','\itset1\Delta','WT MBC',...
    '\itepe1\Delta \rm MBC','\itswi6\Delta \rm MBC',...
    '\itgcn5\Delta \rm MBC','\itset2\Delta \rm MBC',...
    '\itset1\Delta \rm MBC'};
colorAll=[0,0,0; 0,0,1; 1,0,0; .91,.41,.17; .13,.54,.13; .39, .2, .6;...
    0,0,0; 0,0,1; 1,0,0; .91,.41,.17; .13,.54,.13; .39, .2, .6];
numM=length(nameNMBC);
[ points2,faces2,psi,theta] = setup_img( );
p2um=0.16;
f2s=2.5;
[points,faces,edges,neighbors] = TriSphere(3);
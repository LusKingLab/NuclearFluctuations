% Data names and colors.
nameNMBC={'wt','epe1','swi6'};
nameIntercleave=reshape([nameNMBC],1,length(nameNMBC));
nameAll=[nameNMBC];
nameAll2 = {'WT','\it epe1\Delta','\it swi6\Delta'};
colorAll=[0,0,0; 0,0,1; 1,0,0];
numM=length(nameNMBC);
[ points2,faces2,psi,theta] = setup_img( );
p2um=0.16;
f2s=2.5;
[points,faces,edges,neighbors] = TriSphere(3);
namegroup1={'wt','swi6','swi6-sm1'};
namegroup2={'wt_MBC','swi6_MBC','swi6-sm1_MBC'};
nameIntercleave=reshape([namegroup1;namegroup2],1,2*length(namegroup1));
nameAll=[namegroup1,namegroup2];
nameAll2={'WT','\itswi6\Delta','\itswi6-sm1','WT MBC',...
    '\itswi6\Delta \rmMBC','\itswi6-sm1 \rmMBC'};
colorAll=[0, 0, 0;1, 0, 0;1, .41, .71;0, 0, 0;1, 0, 0;1, .41, .71];
% colorAll=[0, 0, 0;1, 0, 0;.25, .88, .82;1, .41, .71];
numM=length(namegroup1);
[ points2,faces2,psi,theta] = setup_img( );
p2um=0.16;
f2s=2.5;
[points,faces,edges,neighbors] = TriSphere(3);
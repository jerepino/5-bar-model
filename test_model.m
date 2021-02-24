q1 = [0.4072    1.5951    0.7321]';
q2 = [1.1393    1.5951         0]';
qda = [0.5; 0.5];
qdda =[1;1];
n=50;
% x = -0.1:0.2/n:0.1;
% y = 0.2:0.1/n:0.3;
p= zeros(n,2);
t0r = p;
t0rd = p;
% figure(5)
% lspb(0.2,0.3,n);
[p(:,1),t0r(:,1),t0rd(:,1)] = lspb(-0.1,0.1,n);
% [p(:,2),t0r(:,2),t0rd(:,2)] = lspb(0.35,0.3,n);
p(:,2) = ones(n,1)*0.3;
for i=1:n
    [q(:,i),~,~] = ikine5(p(i,:)');
    % p_ = fkine5(q(:,i))
    if q(2,i) == 0 || q(2,i) ==pi || q(2,i) ==-pi || q(2,i) ==-2*pi ...
            || q(2,i) ==2*pi || q(5,i)==0 || q(5,i)==pi || q(5,i)==-pi ...
            || q(5,i)==-2*pi || q(5,i)==2*pi
        disp('Type 1 singularity found');
    elseif q(1,i)+q(2,i)-q(4,i)-q(5,i)==0 || q(1,i)+q(2,i)-q(4,i)-q(5,i)==pi ...
            || q(1,i)+q(2,i)-q(4,i)-q(5,i)==-pi || q(1,i)+q(2,i)-q(4,i)-q(5,i)==2*pi ...
            || q(1,i)+q(2,i)-q(4,i)-q(5,i)==-2*pi
        disp('Type 2 singularity found')
    end
%     q(1,i)+q(2,i)-q(4,i)-q(5,i)
    % t0r = ffirstkine5(q,qda);
    [qd(:,i),qda(:,i),~,Ar,B,J,Jinv,Jt,Jta,Jtd,Psit] = ifirstkine5(t0r(i,:)',p(i,:)',q(:,i));
    % qda_

    % t0rd = fsecondkine5(qdda,qd,q,Ar,J);
    [qdd(:,i),qdda(:,i),~,Psit1,Psit2,atp,ad] = isecondkine5(t0rd(i,:)',t0r(i,:)',p(i,:)',q(:,i),qd(:,i),Ar,B,Jinv,Jt,Jta,Jtd,Psit);
    % qdda_
    T(:,i) = idyn5(q(:,i),qd(:,i),qdd(:,i),t0rd(i,:)',J,Jt,Jta,Jtd,Psit);
    qdda_d(:,i) = ddyn5(T(:,i),q(:,i),qd(:,i),qdd(:,i),t0rd(i,:)',J,Jt,Jta,Jtd,Psit,atp,ad);
end
qa = [q(1,:);q(4,:)];

figure(1)
plot(T(1,:),'g--');
hold on;
plot(T(2,:),'k-.');
title('Torque')

figure(2)
plot(qdda_d(1,:),'r-.');
hold on;
plot(qdda(1,:),'b--');
title('Acelerations q11')

figure(3)
plot(qdda(2,:),'m--');
hold on;
plot(qdda_d(2,:),'g-.');
title('Acelerations q21')

figure(4)
plot(q(1,:),'r-.');
hold on;
plot(q(4,:),'b--');

l=0.205;
DH = [  1, 0, l, 0, 0;
        2, 0, l, 0, 0;
        3, 0, 0, 0, 0;];
DH2 =  [4, 0, l, 0, 0;
        5, 0, l, 0, 0;
        6, 0, 0, 0, 0;];
arm1 = SerialLink(DH, 'name', 'arm1');
arm1.base = transl(-0.125, 0, 0);
arm2 = SerialLink(DH2, 'name', 'arm2');
arm2.base = transl(0.125, 0, 0);
% for i=1:4
%     if isnan(q_(1,i))
%         continue;
%     end

figure(5)
plot(p(:,1),p(:,2));
hold on;
arm1.plot(q(1:3,1)')
%trplot(arm1.base * arm1.links(1).A(q(1,1)) )
arm2.plot([q(4:5,1)' 0])
%trplot(arm2.base * arm2.links(1).A(q(4,1)) )

% fkine(arm2,[q(4:5,1)' 0])
% fkine(arm1,q(1:3,1)')
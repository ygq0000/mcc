clc
clear
close all 
f_value = 10000;
load chair3_data.mat; 
run('C:\Users\ygq00\Desktop\3D_cut_250325\vlfeat-0.9.21\toolbox\vl_setup')
YY = chair3(1:15000,:);   %15724 * 3
Q = YY;   QW = YY;
ran = (randperm(size(YY,1)))';  
ZZZ= 4000;
NY = YY(ran(1:ZZZ,:),:) + repmat(rand(1,3) * 0.004 * 2 - 0.004, [ZZZ,1]); 
YY(ran(1:ZZZ,:),:) = NY;
X = YY;
figure(1); plot3(Q(:,1), Q(:,2), Q(:,3), 'ro', X(:,1), X(:,2), X(:,3), 'b*'); axis equal; grid off; axis off; hold off

Q(11000:end,:) = 0;
X(1:6000,:) = 0;
figure(2); plot3(Q(:,1), Q(:,2), Q(:,3), 'ro', X(:,1), X(:,2), X(:,3), 'b*'); axis equal; grid off; axis off; hold off

          x=[-30*pi/180 10*pi/180 40*pi/180 -0.6 1.0 0.5];  
            cx = cos(x(1));cy = cos(x(2));cz = cos(x(3));sx = sin(x(1)); sy = sin(x(2));sz = sin(x(3));
            rr = [cy*cz, sx*sy*cz-cx*sz, cx*sy*cz+sx*sz;
                  cy*sz, cx*cz+sx*sy*sz, cx*sy*sz-sx*cz;
                   -sy,          sx*cy,          cx*cy];
            tt = x(4:6);          
X = X * rr + repmat(tt, [size(X,1), 1]); 
 figure(3); plot3(Q(:,1), Q(:,2), Q(:,3), 'ro', X(:,1), X(:,2), X(:,3), 'b*'); axis equal; grid off; axis off; hold off

 Q = Q'; X = X';
kdtreeQ = vl_kdtreebuild(Q) ;
kdtreeX = vl_kdtreebuild(X) ;
[neighborQ, dQ] = vl_kdtreequery(kdtreeQ, Q, Q, 'numneighbors',5, 'verbose') ;
[neighborX, dX] = vl_kdtreequery(kdtreeX, X, X, 'numneighbors',5, 'verbose') ;
corrpond=find(sum(sort(neighborQ)-sort(neighborX))==0);   %

 QQ = (Q(:,corrpond))';
 XX = (X(:,corrpond))';
figure(4); plot3(QQ(:,1), QQ(:,2), QQ(:,3), 'ro', XX(:,1), XX(:,2), XX(:,3), 'b*'); axis equal; grid off; axis off; hold off

Y2 = QQ; X2 = XX;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Epsilon = 10^-8; 
MaxStep = 30;
LastMinuserror = 10^-6; 
TotalStep = 1;  % last step, if the program is pause, it is the total iterative step
CurrStep = 1;  % current step
Error_of_iter = zeros(MaxStep,1);
Error_of_iter(1) = 10^6;   % The Error
value_t = 0.75;

[nD, ~] = size(X2);   % the number of the model point
[nM, Dim] = size(Y2);  
[nY2] = pcakn(Y2,10);  
nY2 = nY2./sqrt(sum(nY2.^2,2)); 

   R = eye(Dim); %  3*3
   t = (zeros(1,Dim)); % 1×3
   xx = [reshape(R',9,1); t'];
   XT2 = X2 * R + repmat(t, [nM, 1]);   
   XC2 = X2;
  sigma1 = sum(sqrt(sum((XT2 - Y2).^2,2))) * 2/nD;
  sigma2 = sigma1^2;
while (LastMinuserror > Epsilon & TotalStep < MaxStep  & (sigma2 > 0.5));
  
        for i = 1 : nD;
            PPN{i,:} = [  XT2(i,:)  zeros(1,3)  zeros(1,3) 1 0 0; 
                        zeros(1,3)    XT2(i,:)  zeros(1,3) 0 1 0;
                        zeros(1,3)  zeros(1,3)    XT2(i,:) 0 0 1];  
            A(i,:) = nY2(i,:) * PPN{i,:};  
            B(i,:) = nY2(i,:) * (Y2(i,:))';  
            exg(i,:) = ((PPN{i,:} * xx - (Y2(i,:))')' * (nY2(i,:))').^2;
            g(i,:) = exp( -(exg(i,:))./ (2*sigma2) );
        end
           og = find(g>value_t);
           gop = g(og,:);
           xx = (A(og,:) .* repmat(gop,[1,12]))' * A(og,:) \ (A(og,:))' .* repmat(gop',[12,1]) *B(og,:);
           R1 = [(xx(1:3,:))';(xx(4:6,:))'; (xx(7:9,:))'];
           T1 = (xx(10:12,:))';
           R = R1 * R;
           t = (R1 * t' + T1')';
          
          XT2 = (R *  XC2')'; 
          XT2 = XT2 + repmat(t,[nD,1]);     
          
          XT2_X =   (R *  X)' + repmat(t,[size(X,2),1]);    
          figure(8);plot3(Y2(:,1),Y2(:,2),Y2(:,3),'b*', XT2(:,1),XT2(:,2),XT2(:,3),'ro');axis equal;grid on;xlabel('x');ylabel('y');zlabel('z'); axis off;  hold off;
          figure(9);plot3(QW(:,1),QW(:,2),QW(:,3),'b*', XT2_X(:,1),XT2_X(:,2),XT2_X(:,3),'ro');axis equal;grid on;xlabel('x');ylabel('y');zlabel('z'); axis off;  hold off;     
            tempdistance = XT2 - Y2;
            Error_of_iter(CurrStep) = sum(sum(tempdistance.^2, 2))/nD; 
            TotalStep = CurrStep;     
            CurrStep = CurrStep + 1;
            Error_of_last = Error_of_iter(TotalStep);           
                        sigma1 = sum(sqrt(sum((XT2 - Y2).^2,2))) * 2/nD;  
                        sigma2 = sigma1^2;
            if TotalStep == 1
                LastMinuserror = Error_of_iter(TotalStep);
            else
                LastMinuserror = abs(Error_of_iter(TotalStep) - Error_of_iter(TotalStep - 1));
            end         
end
  XT2  = XT2_X;
  figure(10);plot3(Y2(:,1),Y2(:,2),Y2(:,3),'b*', XT2(:,1),XT2(:,2),XT2(:,3),'ro');axis equal;grid off;xlabel('x');ylabel('y');zlabel('z'); axis off;  hold off; 
 
           
           
           
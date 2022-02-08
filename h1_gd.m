% INPUT u0 = noisy image 

% OUTPUT: u(i,j) = denoised image 

% read the noisy data u0 in BMP format

%u0=imread('terminal.bmp');
nb=load('noisybrain.mat')
u0 = nb.Kn
u0=im2double(u0);
[M N]=size(u0);

% visualize the image u0 in Matlab (rescaled) 
imagesc(u0); axis image; axis off; colormap(gray);
%spec_img_log = log(1 + u0);
%imshow(spec_img_log,[]); 


%---------------------------------------------------------------------------
% initialize u by u0 (not necessarily) or by a constant, or by a random image
%---------------------------------------------------------------------------

%lam = 0.01:0.01:0.1
lam = 0.0001:0.0001:0.001

for l = 1:length(lam),
vars = {'u','Energy', 'Fidelity', 'Q'};
clear(vars{:});

u=u0;
[M,N]=size(u);


%----------------------------------------------- 
%           PARAMETERS  
%-----------------------------------------------

% coefficient of the TV norm (needs to be adapted for each image) 
lambda=lam(l); 


%  space discretization 
h=1.; 

%  number of iterations (depends on the image) 
IterMax=500;

% gd rate 
eps=0.0005;


%-----------------------------------------------------------
%     BEGIN ITERATIONS IN ITER 
%-----------------------------------------------------------

l/length(lam)*100

for Iter=1:IterMax, 

	%Iter

    for i=2:M-1,
      for j=2:N-1,

	u(i,j)= (1-4*eps-lambda)*u(i,j) + eps*lambda*u0(i,j) + eps*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1));
    
      end
    end
 
    %------------ FREE BOUNDARY CONDITIONS IN u -------------------
	for i=2:M-1,
          u(i,1)=u(i,2);
          u(i,N)=u(i,N-1);
        end

	for j=2:N-1,
          u(1,j)=u(2,j);
          u(M,j)=u(M-1,j);
        end

        u(1,1)=u(2,2);
        u(1,N)=u(2,N-1); 
        u(M,1)=u(M-1,2);
        u(M,N)=u(M-1,N-1);
%----------------------------------------------------------------------

%---------------- END ITERATIONS IN Iter ------------------------------

%%% Compute the discrete energy at each iteration
en=0.0; q = 0; Fid = 0; 
    for i=2:M-1,
      for j=2:N-1,
      ux=(u(i+1,j)-u(i,j))/h;
      uy=(u(i,j+1)-u(i,j))/h;
      fid=(u0(i,j)-u(i,j))*(u0(i,j)-u(i,j));
      en=en+eps*(ux*ux+uy*uy)+lambda*fid; %
      Fid = Fid + lambda*fid;
      q = q + eps*(ux*ux+uy*uy);
      end
    end
%%% END computation of energy 

Energy(Iter)=en; 
Fidelity(Iter)=Fid;
Q(Iter) = q;

end 


% visualize the image in Matlab (re-scaled)
%figure 
%imagesc(u); axis image; axis off; colormap(gray);
imwrite(u, "h1_gd/brain/h1gd_"+lambda+".jpg", "Quality", 100) 


% use export to save the image in ps or eps format 

% Plot the Energy versus iterations 
%figure 
plot(Energy);legend('Energy/Iterations');
%hold on
%plot(Fidelity);
%hold on
%plot(Q);legend('Energy/Iterations','Fidelity/Iterations','Q/Iterations');
%hold off 
%plot(Energy);legend('Energy/Iterations');
saveas(gcf, "h1_gd/brain/h1gd_Energy_"+lambda+".jpg") 

end
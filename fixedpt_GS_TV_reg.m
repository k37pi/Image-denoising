% Matlab code to denoise an image by the Rudin-Osher-Fatemi model
% Min TV(u) + lambda*(L^2[(f-u)])^2
% lambda = coefficient of the fidelity term 
% lambda needs to be tunned to obtain a sharp denoised image
% h = space discretization step 
% The stationary problem (no time-regularization) is used here
% A fixed point iteration is used 
% The numerical scheme is implicit in the max norm and unconditionally stable
% A stopping criteria can be included, or the "IterMax" must be given. 

% LAST MODIFIED: 2005-04-19 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



%----------------------------------------------- 
%           PARAMETERS  
%-----------------------------------------------

% coefficient of the TV norm (needs to be adapted for each image) 
%lam = 0.5:0.5:5.0
lam = 0.5:0.5:3

for l=1:length(lam),
vars = {'u','Energy', 'Fidelity', 'Q'};
clear(vars{:});

u=u0;
[M,N]=size(u);
    
    
lambda=lam(l); 

%  space discretization 
h=1.; 

%  number of iterations (depends on the image) 
IterMax=100;

%  needed to regularize TV at the origin 
eps=0.0005;


%-----------------------------------------------------------
%     BEGIN ITERATIONS IN ITER 
%-----------------------------------------------------------

l/length(lam)*100

for Iter=1:IterMax, 

	%Iter

    for i=2:M-1,
      for j=2:N-1,

    %-----------computation of coefficients co1,co2,co3,co4---------

	ux=(u(i+1,j)-u(i,j))/h;
	uy=(u(i,j+1)-u(i,j))/h;
        Gradu=sqrt(eps*eps+ux*ux+uy*uy);
        co1=1./Gradu;

        ux=(u(i,j)-u(i-1,j))/h;
	uy=(u(i-1,j+1)-u(i-1,j))/h;
        Gradu=sqrt(eps*eps+ux*ux+uy*uy);
        co2=1./Gradu;

        ux=(u(i+1,j)-u(i,j))/h;
	uy=(u(i,j+1)-u(i,j))/h;
        Gradu=sqrt(eps*eps+ux*ux+uy*uy);
        co3=1./Gradu;

	ux=(u(i+1,j-1)-u(i,j-1))/h;
        uy=(u(i,j)-u(i,j-1))/h;
        Gradu=sqrt(eps*eps+ux*ux+uy*uy);
        co4=1./Gradu;

        co=1.+(1/(2*lambda*h*h))*(co1+co2+co3+co4);

	div=co1*u(i+1,j)+co2*u(i-1,j)+co3*u(i,j+1)+co4*u(i,j-1);

	u(i,j)=(1./co)*(u0(i,j)+(1/(2*lambda*h*h))*div);

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
en=0.0;  
    for i=2:M-1,
      for j=2:N-1,
      ux=(u(i+1,j)-u(i,j))/h;
      uy=(u(i,j+1)-u(i,j))/h;
      fid=(u0(i,j)-u(i,j))*(u0(i,j)-u(i,j));
      en=en+sqrt(eps*eps+ux*ux+uy*uy)+lambda*fid;
      end
    end
%%% END computation of energy 

Energy(Iter)=en; 

end 


% visualize the image in Matlab (re-scaled)
%figure 
%imagesc(u); axis image; axis off; colormap(gray);

imwrite(u, "tv2_gs/terminal/GSTV2_"+lambda+".jpg", "Quality", 100) 

% use export to save the image in ps or eps format 

% Plot the Energy versus iterations 
%figure1 = figure;
%axes1 = axes('Parent',figure1)
%hold(axes1,'all');
plot(Energy);legend('Energy/Iterations');
saveas(gcf, "tv2_gs/terminal/GSTV2_Energy_"+lambda+".jpg") 

end





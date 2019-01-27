%time-optimal problem
%% test

clc 
clear 

%GR1

A=@(t) [0 0;0 0];
B=@(t) [1 0;0 1];
f=@(t) [0;0];

a=0;
b=0;
c=1;
alpha=1/9;
betta=1/25;

r=2;
q1=[56;54];
q2=[74;31];
q3=[93;61];
q4=[41;70];


%GR CROSS

% A=@(t) [0 0;0 0];
% B=@(t) [0 1;3 1];
% f=@(t) [0;0];
% 
% a=0;
% b=4;
% c=3;
% alpha=1/16;
% betta=1;

  %border cross
% r=3;
% q1=[-3;1];
% q2=[5;0];
% q3=[3;12];
% q4=[11;21];

   %X0 in X1
% r=1;
% q1=[0;3];
% q2=[0;-3];
% q3=[3;0];
% q4=[-3;0];

   %X1 in X0
% r=7;
% q1=[0;1];
% q2=[0;-1];
% q3=[1;2];
% q4=[-1;0];

%GR2

% A=@(t) [1 0;2 1];
% B=@(t) [2 12;3 1];
% f=@(t) [-1;10];
% 
% a=0;
% b=0;
% c=1;
% alpha=1;
% betta=1;
% 
% r=2;
% q1=[6;14];
% q2=[4;33];
% q3=[3;22];
% q4=[11;31];


%GR3

% A=@(t) [-3 5;1 1];
% B=@(t) [4 1;3 1];
% f=@(t) [0;10];
% 
% a=-3;
% b=5;
% c=16;
% alpha=1;
% betta=1/4;
% 
% r=2;
% q1=[12;27];
% q2=[17;43];
% q3=[-5;32];
% q4=[26;29];

%GR4

% A=@(t) [0 -1;1 0];
% B=@(t) [0 1;0 0];
% f=@(t) [0;0];
% 
% a=0;
% b=0;
% c=1;
% alpha=1;
% betta=1;
% 
% r=2;
% q1=[6;3];
% q2=[10;3];
% q3=[8;7];
% q4=[8;-1];


%GR5

% A=@(t) [0 1;-2 3];
% B=@(t) [0 -cos(t^2)/2;sin(t) 0];
% f=@(t) [0.1*t;0];
% 
% a=0;
% b=0;
% c=1;
% alpha=1;
% betta=1/4;
% 
% r=1;
% 
% q1=[-5;-2];
% q2=[-7;2];
% q3=[-2.5;1];
% q4=[-5;0];


%GR6
% A=@(t) [15 10;-10 15];
% B=@(t) [0 1;sin(t/10) 0];
% f=@(t) [0;30*cos(30*t)];
% 
% a=0;
% b=0;
% c=16;
% alpha=1;
% betta=1;
% 
% r=1;
% 
% q1=[8;13];
% q2=[10;11];
% q3=[11;3];
% q4=[7;2];
% 

q=[q1,q2,q3,q4];

%% main part
% STREIGHT TASK
T=3;
tstep=0.01;
N=25;

cont=1;
qc=convQ(q1,q2,q3,q4);

%CROSS?

i=1;
angle=linspace(0,2*pi,N);
convX1= @(L) max([dot(L,q(:,1)),dot(L,q(:,2)),dot(L,q(:,3)),dot(L,q(:,4))]);

cross=1;
while((i<=N)&&(cross))
  psi=[cos(angle(i)),sin(angle(i))];
  if(convX1(psi)<0)
    cross=0;  
  end    
  i=i+1;
end    

if((qc(1,1)^2+qc(2,1)^2<=r^2)||(cross))
  cross=1;
else
  cross=0;  
end    


k=1;
while((k<length(qc(1,:))-1)&&(~cross))
  a1=qc(1,k);
  a2=qc(2,k);
  b1=qc(1,k+1);
  b2=qc(2,k+1);
  leftX=min(a1,b1);
  rightX=max(a1,b1);
  leftY=min(a2,b2);
  rightY=max(a2,b2);
  if((((a1-b1)^2+(a2-b2)^2)*r^2-(a1*b2-a2*b1)^2)>=0)  
    yc=((b1-a1)*(a2*b1-a1*b2)+(a2-b2)*(((a1-b1)^2+(a2-b2)^2)*r^2-(a1*b2-a2*b1)^2)^(1/2))/((a1-b1)^2+(a2-b2)^2);          
    xc=((b2-a2)*(a1*b2-a2*b1)+(a1-b1)*(((a1-b1)^2+(a2-b2)^2)*r^2-(a1*b2-a2*b1)^2)^(1/2))/((a1-b1)^2+(a2-b2)^2);
    if((xc>=leftX)&&(xc<=rightX)&&(yc>=leftY)&&(yc<=rightY))
      cross=1;
    end    
    xc=((b2-a2)*(a1*b2-a2*b1)-(a1-b1)*(((a1-b1)^2+(a2-b2)^2)*r^2-(a1*b2-a2*b1)^2)^(1/2))/((a1-b1)^2+(a2-b2)^2);
    yc=((b1-a1)*(a2*b1-a1*b2)-(a2-b2)*(((a1-b1)^2+(a2-b2)^2)*r^2-(a1*b2-a2*b1)^2)^(1/2))/((a1-b1)^2+(a2-b2)^2); 
    if((xc>=leftX)&&(xc<=rightX)&&(yc>=leftY)&&(yc<=rightY))
      cross=1;
    end   
  end
  k=k+1;
 end

if(cross)
  ts=fill(qc(1,:),qc(2,:),'c');
  hold on
  angle=linspace(0,2*pi,N);
  is=plot(r*cos(angle),r*sin(angle),'k');
  hold off
  legend([is,ts],'initial set','target set');
  title('Initial and target sets have common points. Time-optimal problem: T = 0');
  xlabel('x1');
  ylabel('x2');
else  
  while cont
    tstepg=tstep;
    Tg=T;  
    figM=figure('Name','Main plot'); 
    fmain=subplot(1,1,1);  
    
    ts=fill(qc(1,:),qc(2,:),'c');
    hold(fmain,'on');
    
    optU=figure('Name','Optimal control');
    ou=subplot(1,1,1);
    
    supP=@(psi) [a+psi(1)*(betta*c/(alpha*(betta*(psi(1)^2)+alpha*(psi(2)^2))))^(1/2);
                 b+psi(2)*(alpha*c/(betta*(betta*(psi(1)^2)+alpha*(psi(2)^2))))^(1/2)];
    
    phi=linspace(0,2*pi,2*N);
    psi=[cos(phi);sin(phi)];
    hold(ou,'on')
    xlabel(ou,'u1');
    ylabel(ou,'u2');
    title(ou,'Optimal control');
    axis(ou,[a-(c/alpha)^(1/2)-0.3 a+(c/alpha)^(1/2)+0.3 b-(c/betta)^(1/2)-0.3 b+(c/betta)^(1/2)+1.5]);
    precontrol=supP(psi(:,end));
    for i=1:2*N
      control=supP(psi(:,i)); 
      pCont=plot(ou,[precontrol(1) control(1)],[precontrol(2) control(2)],'b');
      precontrol=control;
    end
    
    fig=figure('Name','Trajectories');
    f1=subplot(3,2,1);
    f2=subplot(3,2,2);
    f3=subplot(3,2,3);
    f4=subplot(3,2,4);
    f5=subplot(3,2,5);
    f6=subplot(3,2,6);
    hold(f1,'on');
    hold(f2,'on');
    hold(f3,'on');
    hold(f4,'on');
    hold(f5,'on');
    hold(f6,'on');
    
    mintime=-1;
    angle=linspace(0,2*pi,N);
    turn=@(psi) [cos(psi) -sin(psi); sin(psi) cos(psi)];%turn on psi
    is=plot(fmain,r*cos(angle),r*sin(angle),'k');
    for i=1:length(angle)
      t=0;
      psi0=[cos(angle(i));sin(angle(i))];
      x0=r*psi0;
      final=[];
      u0=supP((B(t))'*psi0);
      while((isempty(final))&&(t<T))
        H=@(t,y) systemH(A,B,f,u0,t,y);
        options=odeset('Events',@(t,y) AttainableSet(y,q));
        [t1,y,te,ye,ie]=ode45(H,[t,t+tstep],[x0;psi0],options);
        x=[y(end,1);y(end,2)];
        final=te;
      
        psi=[y(end,3);y(end,4)];
        u=supP((B(t1(end)))'*psi);
        
        plot(f1,[t t1(end)],[x0(1),x(1)],'c');
        plot(f2,[t t1(end)],[x0(2),x(2)],'c');
        plot(f3,[t t1(end)],[u0(1),u(1)],'c')
        plot(f4,[t t1(end)],[u0(2),u(2)],'c');
        plot(f5,[t t1(end)],[psi0(1),psi(1)],'c');
        plot(f6,[t t1(end)],[psi0(2),psi(2)],'c');

        t=t1(end);
        tr=plot(fmain,[x0(1),x(1)],[x0(2),x(2)],'b');
        x0=x;
        u0=u;
        psi0=psi;
      end
      if((t<mintime)||(mintime==-1))
        mintime=t;
        minnumb=i;
      end    
    end
  
% NOT OPTIM READY
    if(mintime>=T)
      notoptim=1;
      disp('The problem is unsolvable');
      legend(fmain,[is,ts,tr],'initial set','target set','trajectory');
      txt='Time-optimal problem has no solution';
      title(fmain,txt);
      xlabel(fmain,'x1');
      ylabel(fmain,'x2');
    else 
      notoptim=0;  
      phi=angle(minnumb);
      psi=[cos(phi);sin(phi)];
      x0=r*psi;
      final=[];
      t=0;
      mintrack=zeros(6,ceil((T)/tstep));
      k=1;
      u=supP((B(t))'*psi);
      mintrack(:,k)=[x0;psi;u];
      while(isempty(final))
        H=@(t,y) systemH(A,B,f,u,t,y);
        options=odeset('Events',@(t,y) AttainableSet(y,q));
        [t1,y,te,ye,ie]=ode45(H,[t,t+tstep],[x0;psi],options);
        x=[y(end,1);y(end,2)];
        final=te;
        t=t1(end);
      
        psiNext=[y(end,3);y(end,4)];
      
        pVal=plot(ou,u(1),u(2),'or');
       
        u=supP((B(t1(end)))'*psiNext);
        k=k+1;
        mintrack(:,k)=[x;psiNext;u];
      
        opttr=plot(fmain,[x0(1),x(1)],[x0(2),x(2)],'r');
        if(isempty(final))
          x0=x;
          psi=psiNext;
        end
      end
      mintrack=mintrack(:,1:k);
      legend(ou,[pCont,pVal],'Boundary of the optimal control set','Real values of optimal control'); 
    %hold(ou,'off');
    
%Error of transversality conditions
      if(length(qc(1,:))==1)
        eps=0;  
      else   
        if(length(qc(1,:))>=3)  
          eps=0.001;
          [~,j]=min([dot(psi,q(:,1)),dot(psi,q(:,2)),dot(psi,q(:,3)),dot(psi,q(:,4))]);
          i=find((qc(1,:)-q(1,j)+qc(2,:)-q(2,j))<=eps,1);
        
         q1=[qc(:,end-1) qc];
         qnorm=[q1(:,i) q1(:,i+1) q1(:,i+2)];
        else  
          qnorm=qc;  
        end  
        x0norm=x0;
        xnorm=x;
        eps=drawNorm(fmain,psi,qnorm,x0norm,xnorm,1);%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!I'M HERE!!!!!!!
        %convX1= @(L) max([dot(L,q(:,1)),dot(L,q(:,2)),dot(L,q(:,3)),dot(L,q(:,4))]);
      end  
      str=['Error of transversality conditions = ',num2str(eps)];
      disp(str);
%Error end
      
%DRAW
      legend(fmain,[is,ts,tr,opttr],'initial set','target set','trajectory','optimal trajectory');
      txt=['Time-optimal problem: T= ',num2str(mintime,3)];
      disp(txt);
      title(fmain,txt);
      xlabel(fmain,'x1');
      ylabel(fmain,'x2');
      
      Mx=max([r,q(1,:)]);
      mx=min([-r,q(1,:)]);
      My=max([r,q(2,:)]);
      my=min([-r,q(2,:)]);
      ots=3;
      if(Mx-mx>My-my)
        h=Mx-mx;  
        center=(My-my)/2;
        axis(fmain,[mx-ots Mx+ots center-h/2-ots center+h/2+ots]);
      else    
        h=My-my;
        center=(Mx-mx)/2;
        axis(fmain,[center-h/2-ots center+h/2+ots my-ots My+ots]);
      end    

      if(~notoptim)      
%OPTIMAL TRAJECTORY
        tvec=0:tstep:tstep*(length(mintrack(1,:))-1);
        plot(f1,tvec,mintrack(1,:),'r');
        title(f1,'X-component of optimal trajectory');
        xlabel(f1,'t');
        ylabel(f1,'x1');
     
%     subplot(2,1,2);
        plot(f2,tvec,mintrack(2,:),'r');
        title(f2,'Y-component of optimal trajectory');
        xlabel(f2,'t');
        ylabel(f2,'x2');
%OPTIMAL CONTROL
%     figure
%     subplot(2,1,1);
        plot(f3,tvec,mintrack(5,:),'r');
        title(f3,'X-component of optimal control');
        xlabel(f3,'t');
        ylabel(f3,'u1');
%     
%     subplot(2,1,2);
        plot(f4,tvec,mintrack(6,:),'r');
        title(f4,'Y-component of optimal control');
        xlabel(f4,'t');
        ylabel(f4,'u2');
%ADJOINT VARIABLES
%     figure
%     subplot(2,1,1);
        plot(f5,tvec,mintrack(3,:),'r');
        title(f5,'X-component of adjoint variable');
        xlabel(f5,'t');
        ylabel(f5,'psi1');
    
%     subplot(2,1,2);
        plot(f6,tvec,mintrack(4,:),'r');
        title(f6,'Y-component of adjoint variable');
        xlabel(f6,'t');
        ylabel(f6,'psi2');
      end
      axis([f1,f2,f3,f4,f5,f6],[0 T -inf inf]);
  
      hold(f1,'off');
      hold(f2,'off');
      hold(f3,'off');
      hold(f4,'off');
      hold(f5,'off');
      hold(f6,'off');
  
      
%LOCAL OPTIMIZATION     
      angle=[angle(end) angle angle(1)];
      optLoc=1;
      NLoc=1;
      tstepLoc=tstep;
      TLoc=T;
      phi=angle(minnumb+1);
      while optLoc
        tryagain=1;
        tryes=1;
        %N=NLoc;
        tstep=tstepLoc;
        T=TLoc;
        while((tryagain)&&(tryes<=3))
          promt='Do you want to change parametres for local optimization? y/n? ';
          answ=input(promt);
          y='y';
          n='n';
          if(answ==y)
            promt=['Current value of time step is ',num2str(tstepLoc),'; Enter new value: '];
            tstepLoc=input(promt);
            promt=['Current number of points between the closest points of global optimisation is ',num2str(NLoc),'; Enter new value:'];
            NLoc=input(promt)+2;
            promt=['Current time limit is ',num2str(TLoc),'; Enter new value:'];
            TLoc=input(promt);
            tryagain=0;
          else
            if(answ~=n)
              disp('answer should be ''y'' or ''n''');
              tryagain=1;
            else
              tryagain=0;  
            end    
            optLoc=0;    
          end
          tryes=tryes+1;
        end 
      
        if(optLoc)
%HIDE PREVIOUS OPTIMAL TRAJECTORY          
          psi=[cos(phi);sin(phi)];
          x0=r*psi;
          final=[];
          t=0;
          while(isempty(final))
            u=supP((B(t))'*psi);
            H=@(t,y) systemH(A,B,f,u,t,y);
            options=odeset('Events',@(t,y) AttainableSet(y,q));
            [t1,y,te,ye,ie]=ode45(H,[t,t+tstep],[x0;psi],options);
            x=[y(end,1);y(end,2)];
            final=te;
            t=t1(end);
            opttr=plot(fmain,[x0(1),x(1)],[x0(2),x(2)],'m');
            x0=x;
            if(isempty(final))
              psi=[y(end,3);y(end,4)];
            end  
          end
          if(length(qc(1,:))>1)
            drawNorm(fmain,psi,qnorm,x0norm,xnorm,0);
          end  
          
%END HIDE 
      
          angles=linspace(angle(minnumb),angle(minnumb+2),NLoc);
          premintime=mintime;
          mintime=-1;
          for i=1:length(angles)
            t=0;
            angeLoc=angles(i);
            psi=[cos(angeLoc);sin(angeLoc)];
            x0=r*psi;
            final=[];
            u=supP((B(t))'*psi);
            xtrack=zeros(6,ceil((TLoc)/tstepLoc));
            k=1;
            xtrack(:,k)=[x0;psi;u];
            while((isempty(final))&&(t<TLoc))
              H=@(t,y) systemH(A,B,f,u,t,y);
              options=odeset('Events',@(t,y) AttainableSet(y,q));
              [t1,y,te,ye,ie]=ode45(H,[t,t+tstepLoc],[x0;psi],options);
              psi=[y(end,3);y(end,4)];
              u=supP((B(t1(end)))'*psi);
              x=[y(end,1);y(end,2);y(end,3);y(end,4);u];
              final=te;
              t=t1(end);
              k=k+1;
              xtrack(:,k)=x;
              x0=[y(end,1);y(end,2)];
            end
            xtrack=xtrack(:,1:k);
            if((t<mintime)||(mintime==-1))
              mintime=t;
              mintrack=xtrack;
              phi=angeLoc;
            end    
          end
          opttrnew=plot(fmain,mintrack(1,:),mintrack(2,:),'r');
          
          if(mintime<TLoc)
            if(length(qc(1,:))==1)
              eps=0;  
            else    
              psi=[mintrack(3,end-1);mintrack(4,end-1)];
              if(length(qc(1,:))>=3)
                eps=0.001;
                [~,j]=min([dot(psi,q(:,1)),dot(psi,q(:,2)),dot(psi,q(:,3)),dot(psi,q(:,4))]);
                i=find((qc(1,:)-q(1,j)+qc(2,:)-q(2,j))<=eps,1);
       
                q1=[qc(:,end-1) qc];
                qnorm=[q1(:,i) q1(:,i+1) q1(:,i+2)];
              else
                qnorm=qc;  
              end  
              x0norm=[mintrack(1,end-1);mintrack(2,end-1)];
              xnorm=[mintrack(1,end);mintrack(2,end)];
              eps=drawNorm(fmain,psi,qnorm,x0norm,xnorm,1);
            end  
            str=['Error of transversality conditions = ',num2str(eps)];
            disp(str);
          end
          
          
          legend(fmain,[is,ts,tr,opttr,opttrnew],'initial set','target set','trajectory','previous optimal trajectory','new optimal trajectory');
          txt=['Time-optimal problem: T= ',num2str(mintime,3)];
          title(fmain,txt);
          disp(txt);
          %hold(fmain,'off');
        end    
      end
    end 
  
%CHANGE PARAMETRES
    %hold(fmain,'off');
    tryagain=1;
    tryes=1;
    while((tryagain)&&(tryes<=3))
      promt='Do you want to evaluate this task with other global paramaters of accuracy? y/n? ';
      answ=input(promt);
      y='y';
      n='n';
      if(answ==y)
        cont=1;
        promt=['Current value of time step is ',num2str(tstepg),'; Enter new value: '];
        tstep=input(promt);
        promt=['Current number of points of geometric split is ',num2str(N),'; Enter new value:'];
        N=input(promt);
        promt=['Current time limit is ',num2str(Tg),'; Enter new value:'];
        T=input(promt);
        tryagain=0;
      else
        if(answ~=n)
          disp('answer should be ''y'' or ''n''');
          tryagain=1;
        else
          tryagain=0;  
        end    
        cont=0;    
      end
      tryes=tryes+1;
    end  
  end
end 
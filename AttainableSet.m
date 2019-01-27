function [value,isterminal,direction]=AttainableSet(y,q)

  convX1= @(L) max([dot(L,q(:,1)),dot(L,q(:,2)),dot(L,q(:,3)),dot(L,q(:,4))]);
  
  eps=-0.1;
  value=0;
  angle=linspace(0,pi,20);
  i=1;
  while((i~=length(angle))&&(value~=1))
    psi=[cos(angle(i)),sin(angle(i))];
    if(((dot(psi,y(1:2))-eps>convX1(psi)))||((dot(psi,y(1:2)))+eps<-convX1(-psi)))
      value=1;
    end
    i=i+1;
  end 
  isterminal=1;
  direction=0;
end
function [q]=convQ(q1,q2,q3,q4)
  N=90;
  angle=linspace(0,2*pi,N);
  points=[q1,q2,q3,q4];
  q=[];
  for i=1:N
    psi=[cos(angle(i)),sin(angle(i))];  
    [~,ind]=max([dot(psi,q1),dot(psi,q2),dot(psi,q3),dot(psi,q4)]);
    if isempty(q)
      q=points(:,ind);
    else
      if(points(:,ind)~=q(:,end))
        q=[q,points(:,ind)];
      end  
    end  
  end    

  
  
%   q=[q1,q2,q3,q4];
%   k=@(x1,x2) (x2(2)-x1(2))/(x2(1)-x1(1));
%   b=@(x1,x2) x2(2)-k(x1,x2)*x2(1);
%   p1=@(x1,x2,z) (z(1)+k(x1,x2)*(z(2)-b(x1,x2)))/(k(x1,x2)^2+1);
%   p2=@(x1,x2,z) (z(2)*(k(x1,x2))^2+z(1)*k(x1,x2)+b(x1,x2))/(k(x1,x2)^2+1);
%   distlp=@(ln1,ln2,pnt) ((pnt(1)-p1(ln1,ln2,pnt))^2+(pnt(2)-p2(ln1,ln2,pnt))^2)^(1/2);
%   if((q1==q2)&&(q2==q3)&&(q3==q4))%POINT
%     q=q1;
%   else
%     if((distlp(q1,q2,q3)==0)&&(distlp(q1,q2,q4)==0))%LINE
%       [~,left]=min([norm2(q1),norm2(q2),norm2(q3),norm2(q4)]);
%       [~,right]=max([norm2(q1),norm2(q2),norm2(q3),norm2(q4)]);
%       q=[q(left),q(right)];
%     else
%       if((distlp(q1,q2,q3)==0)||(distlp(q1,q2,q4)==0)||(distlp(q2,q3,q4)==0))%TRIANGLE
%         if(distlp(q1,q2,q3)==0)
%           [~,left]=min([norm2(q1),norm2(q2),norm2(q3)]);
%           [~,right]=max([norm2(q1),norm2(q2),norm2(q3)]);
%           q=[q(left),q(right),q4];   
%         end  
%         if(distlp(q1,q2,q4)==0)
%           q=[q1,q2,q4];  
%           [~,left]=min([norm2(q1),norm2(q2),norm2(q4)]);
%           [~,right]=max([norm2(q1),norm2(q2),norm2(q4)]);
%           q=[q(left),q(right),q3];   
%         end 
%         if(distlp(q2,q3,q4)==0)
%           q=[q2,q3,q4];  
%           [~,left]=min([norm2(q2),norm2(q3),norm2(q4)]);
%           [~,right]=max([norm2(q2),norm2(q3),norm2(q4)]);
%           q=[q(left),q(right),q1];   
%         end 
%       else %QUADRANGLE
%         if()  
%       end    
%     end  
%   end    
end 

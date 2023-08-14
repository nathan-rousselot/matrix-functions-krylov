function [readyPoly , polyn , ar] = hermitePoly( ar2 )
%HERMITEPOLY return the polynom calculated by Hermite method.
%   input: ar2 is a matrix. The 1st row is the x values, the 2nd row is the
%   f(x) values. The 3th is the first derivative, f'(x), if exist. the 4th
%   row is the second derivative etc.
%   If any derivative for a point does not exist, e.g. f'(x4), but the row
%   exist because f'(x!=x4) is known, then the derivative value in the
%   matrix MUST be NaN.
%   example: x1=1 x2=2 x3=3 f(x1)=10 f(x2)=20 f(x3)=30 f'(x1)=100 f'(x3)=300
%   the input in this case would be [1 2 3; 10 20 30; 10 NaN 300]
%
%   output: readyPoly - a vector of interpolation polynomial redy for use
%   in polyval function. [An An-1 .... A0]
%           polyn - the parameters that taken from the table before
%           multiplying An(x-x0)(x-x1)...(x-xn-1). can be used to check
%           parameters if calculated by hand.
%           ar - table divided differences.
%   version 1.2
%   wroten by Or Werner, BGU, Israel. For cantact please use mathworks file
%   exchange.

field1='val';
field2='xmin';
field3='xmax';
s=struct(field1,NaN,field2,[],field3,[]);
n=length(ar2(1,:));
levels=length(ar2(:,1))-1;

for ii=1:levels
ar{ii}=s;
end

for ii=1:n %running on xi->xn
    
    xi=ar2(1,ii); %the current x.
    
    for jj=1:levels %running on f , f' , f''
        if ~(isnan(ar2(jj+1,ii))) % if the level (f, f'..) is having a value, building another same field.
            for kk=1:jj %the two "for"s are buildind a piramys style.
                ar{kk}(end+1)=struct(field1,ar2(kk+1,ii),field2,xi,field3,xi);     

            end
        end
    end
    if ii<n
      for jj=2:levels %fill the empty places
          ar{jj}(end+1:length(ar{1}))=struct(field1,NaN,field2,[],field3,[]);
      end
    elseif ii==n
     for jj=2:levels
         ar{jj}(end+1:(length(ar{1})-jj+1))=struct(field1,NaN,field2,[],field3,[]);
     end
    end

end
for ii=1:length(ar)
    ar3{ii}=ar{ii}(2:end); %arranges the varaibles.
end
ar=ar3; %moving the arrays to the varaible ar.

for ii=2:levels
    for jj=1:(length(ar{ii-1})-1)
       if isnan(ar{ii}(jj).val) %calculating the derivatives where they are not given.
           ar{ii}(jj).val=(ar{ii-1}(jj).val-ar{ii-1}(jj+1).val)./(ar{ii-1}(jj).xmin-ar{ii-1}(jj+1).xmax);
           ar{ii}(jj).xmin=ar{ii-1}(jj).xmin;
           ar{ii}(jj).xmax=ar{ii-1}(jj+1).xmax;
       end
    end
end

while length(ar{end})>1 %the loop cuts the un-needed end of every array.
    lengthOfCol=length(ar{end});
    lastCol=length(ar);
    for ii=1:(lengthOfCol-1)
        ar{lastCol+1}(ii)=struct(field1,((ar{lastCol}(ii).val-ar{lastCol}(ii+1).val)./(ar{lastCol}(ii).xmin-ar{lastCol}(ii+1).xmax)),field2,ar{lastCol}(ii).xmin,field3,ar{lastCol}(ii+1).xmax);
    end    
end
polyn=[];
for ii=1:length(ar)
    polyn=[ar{ii}(1).val,polyn]; %creating the array of the polynomial's parameters as taken from the table like calculating by hand.
end


%calculating a polynomial ready for use, without need of (x-x0)(x-x1)
%etc...

%The main idea is to make convolution of (x-x0)...(x-xn) step by step to
%create array of polynomials. Than make all the polynomails to be the
%same range ([0 0 0 1 1] instead of [1 1]). And than multiply every polynomial
%with it's parameter (a0, a1...) and sum. that will be a polynomial ready
%to use with polyval.

conpoly=cell(1,length(polyn));
readyPoly=0;
flippoly=fliplr(polyn);
for ii=1:length(polyn)
    conpoly{ii}=1;
   for  jj=1:(ii-1)
       conpoly{ii}=conv(conpoly{ii},[1,-ar{1}(jj).xmin]); %multiplying (x-x0)(x-x1)... etc.
   end
   conpoly{ii}=flippoly(ii)*[zeros(1,(length(polyn)-length(conpoly{ii}))),conpoly{ii}]; %Arrange all the polynomials to be of the same length. Multyplying with parameter from the table.
   readyPoly=readyPoly+conpoly{ii}; %sum the sub-polynomials.
end
end
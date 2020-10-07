function c = minus(a,b)

   % Convert possible non-PhysicalData "a" to a Value
   if ~isa(a,'PhysicalData')
      try
         a = newvalue(a);          
      catch ME
         % Error if b does not match the expected data format of newvalue.
         rethrow(ME);
      end
   end

   % Convert possible non-PhysicalData "b" to a Value
   if ~isa(b,'PhysicalData')
      try
         b = newvalue(b);          
      catch ME
         % Error if b does not match the expected data format of newvalue.
         rethrow(ME);
      end
   end
   

   % Dispatch the minus operation to the most specialized class
   if strcmp(class(b),'Value') % || more general object
      % Call a - b method of "a"
      c = lhs_minus(a,b);
   else % "b" is a more specialized object
      % Call a - b method of "b"
      c = rhs_minus(b,a);
   end
   
end
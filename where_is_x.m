 function where_is_x = where_is_x( x, x_mesh )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  find which interval in x_mesh into which x(j) falls.
%  on entry, x is a row vector.
%  on exit, where_is_x is a row vector ( length(x) )
%  with entries giving the interval number in x_mesh
% ---------------------------------------------------



%  clunky but simple search
%  -------------------------
       N_x        = length( x );

       where_is_x = NaN( N_x, 1 );
       
      for j=1:N_x

   %    index = find( x_mesh <= x(j) );
   %     where_is_x(j) = index(end)

         where_is_x(j) = find(x_mesh <= x(j), 1, 'last');

      end  %  for j=1:N_x ...
      
      where_is_x = where_is_x';
%
%
%


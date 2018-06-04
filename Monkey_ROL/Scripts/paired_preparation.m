
% load 1000x2 matrix

indexes = [];

for i = 1:length(onsets)
    
   if ~isnan(onsets(i,1)) && ~isnan(onsets(i,2))
       
       indexes = [indexes, i];

       
   end
    
end

indexes = indexes';

new_onsets = onsets(indexes, :);

a = new_onsets(:,1);
b = new_onsets(:,2);
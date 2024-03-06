% Load the JSON data from a file
jsonFileName = 'data.interpolate.json'; % Replace with the actual JSON file path
jsonText = fileread(jsonFileName);

% Parse the JSON data using jsondecode
jsonData = jsondecode(jsonText);

tensorData = tensor(jsonData);


for rank = 2:10
   M = gcp_opt(tensorData,rank,'type','beta (0.5)');
   csvwrite(append("GCP_R", int2str(rank), ".csv"),M.U{1});
end






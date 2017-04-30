folder_path = '../data';
new_folder_path = '../html/pics';
mkdir(new_folder_path);
test_data = dir(fullfile(folder_path,'*.bmp'));

for i=1:size(test_data, 1)
    file_name = test_data(i).name;
    file_path = fullfile(folder_path, file_name);
    new_file_name = strrep(file_name, '.bmp', '.jpg');
    new_file_path = fullfile(new_folder_path, new_file_name);
    im = imread(file_path);
    imwrite(im, new_file_path); 
end
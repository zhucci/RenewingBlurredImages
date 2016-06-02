function varargout = deconv(varargin)
%DECONV M-file for deconv.fig
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @deconv_OpeningFcn, ...
                   'gui_OutputFcn',  @deconv_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before deconv is made visible.
function deconv_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.

% Choose default command line output for deconv
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = deconv_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu_fcn_im.
function popupmenu_fcn_im_Callback(hObject, eventdata, handles)
imageName = get(handles.popupmenu_fcn_im,'String');
PopId = get(handles.popupmenu_fcn_im, 'Value');

J =  imread(imageName{PopId});
handles.cur_normal_image=J;
handles.cur_blur_image = J;
axes(handles.axes_plot);
imshow(handles.cur_normal_image);
  switch get(handles.text3,'Visible')
   
       case 'off'
           set(handles.text3,'Visible','on');
           set(handles.popupmenu_blur,'Visible','on');
           set(handles.text4,'Visible','on');
           set(handles.text6,'Visible','on');
           set(handles.edt_Disk,'Visible','on');
           set(handles.edt_noise,'Visible','on');
           set(handles.btn_Conv, 'Visible','on');
            set(handles.text15, 'Visible','on');
      case 'on'
           set(handles.show_FFT,'Visible','off');
           
  end
  set(handles.call_back,'String','Image has been opened!');
guidata(handles.figure_graphic, handles);

%set(handles.axes_plot,'UIContextMenu','CMenu');
 
function popupmenu_fcn_im_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
end


% --- Executes on selection change in popupmenu_blur.
function popupmenu_blur_Callback(hObject, eventdata, handles)
method_name = get(handles.popupmenu_blur,'String');
PopId = get(handles.popupmenu_blur, 'Value');
switch method_name{PopId}
    case 'Blur'
        set(handles.text4,'Visible','on');
        set(handles.text6,'Visible','on');
        set(handles.edt_Disk,'Visible','on');
        set(handles.edt_noise,'Visible','on');
        set(handles.btn_Conv, 'Visible','on');
        set(handles.show_FFT,'Visible','off');
    case 'No Blur'
      
         handles.cur_blur_image = handles.cur_normal_image;
        set(handles.text4,'Visible','off');
        set(handles.text6,'Visible','off');
        set(handles.edt_Disk,'Visible','off');
        set(handles.edt_noise,'Visible','off');
        set(handles.btn_Conv, 'Visible','off');
        set(handles.btn_save_degrad, 'Visible','off');
         set(handles.text15, 'Visible','off');
         
          if(strcmp(get(handles.btn_renew,'Visible'),'off'))
         set(handles.txt_rest,'Visible','on');
        set(handles.text7,'Visible','on');
        set(handles.popupmenu_restor,'Visible','on');
        set(handles.button_FFT_renew,'Visible','on');
        set(handles.text_noise_adv,'Visible','on');
        set(handles.edt_rad,'Visible','on');
        set(handles.edt_nsr,'Visible','on');
        set(handles.btn_renew,'Visible','on');
        set(handles.checkbox_fast,'Visible','on');
       end
         guidata(gcbo,handles);
end
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_blur contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_blur


% --- Executes during object creation, after setting all properties.
function popupmenu_blur_CreateFcn(hObject, eventdata, handles)
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sl_noise_Callback(hObject, eventdata, handles)
% hObject    handle to sl_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val =get(hObject,'Value');
btn_renew_Callback(hObject, eventdata, handles, val);

function sl_noise_CreateFcn(hObject, eventdata, handles)
% Hint: slider controls usually have a light gray axes_background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function open_new_image(hObject, eventdata, handles)
[fname, pname] = uigetfile('*.png');
if fname ~=0
    n = numel(get(handles.popupmenu_fcn_im, 'String'));
    mass = [];
    mass = get(handles.popupmenu_fcn_im, 'String');
   fullname = strcat(pname, fname);
   mass{end+1} = fullname;
   set(handles.popupmenu_fcn_im, 'String', mass);
   set(handles.popupmenu_fcn_im, 'Value',n+1);
   popupmenu_fcn_im_Callback(hObject, eventdata, handles);
    set(handles.text12, 'Visible','off');
   switch get(handles.text3,'Visible')
       case 'off'
           set(handles.text3,'Visible','on');
           set(handles.popupmenu_blur,'Visible','on');
           set(handles.text4,'Visible','on');
           set(handles.text6,'Visible','on');
           set(handles.edt_Disk,'Visible','on');
           set(handles.edt_noise,'Visible','on');
           set(handles.btn_Conv, 'Visible','on');
          set(handles.text15, 'Visible','on');
        
           
   end
else  
    set(handles.call_back,'String','Image not found(Error 003)');
end

% --- Executes on button press in btn_show_im.
function btn_show_im_Callback(hObject, eventdata, handles)
    open_new_image(hObject, eventdata, handles);
 

% --- Executes on button press in btn_Conv.
function btn_Conv_Callback(hObject, eventdata, handles)
method_name = get(handles.popupmenu_blur,'String');
PopId = get(handles.popupmenu_blur, 'Value');

switch method_name{PopId}
    case 'Blur'
    disk = str2num(get(handles.edt_Disk, 'String'));

    noise_var = (str2num(get(handles.edt_noise, 'String'))*power(10,-6));
    noise_mean = 0;
    m = mod(disk, 1);
    if disk<=0 
    set(handles.call_back,'String','Radius illegal(Error 003)');
    return;
    elseif isempty(disk)
         disk = 10;
         set(handles.edt_Disk,'String',num2str(disk));   
    elseif m ~=0
        set(handles.call_back,'String','Radius must be integer(Error 003)');
    return;
    end
     if isempty(noise_var)
         noise_var = 0.0;
        set(handles.edt_noise,'String',num2str(0));
     elseif noise_var<0
        set(handles.call_back,'String','Noise illegal(Error 003)');
        return;
    end
    PSF = fspecial('disk',disk);
    try
    Blurred = imfilter(handles.cur_normal_image, PSF, 'circula','conv');
    handles.cur_blur_image= imnoise(Blurred, 'gaussian' , noise_mean, noise_var);
    catch
    set(handles.call_back,'String','Out of memory(Error 001)');
    return;
    end;
    axes(handles.axes_plot);
    imshow(handles.cur_blur_image);
     set(handles.call_back,'String','Image has been blurred');
    guidata(handles.figure_graphic, handles);
       
end
        set(handles.btn_save_degrad,'Visible','on');
     if(strcmp(get(handles.btn_renew,'Visible'),'off'))
        set(handles.text7,'Visible','on');
        set(handles.popupmenu_restor,'Visible','on');
        set(handles.button_FFT_renew,'visible','on');
        set(handles.text_noise_adv,'Visible','on');
        set(handles.txt_rest,'Visible','on');
        set(handles.edt_rad,'Visible','on');
        set(handles.edt_nsr,'Visible','on');
        set(handles.btn_renew,'Visible','on');  
        set(handles.checkbox_fast,'Visible','on');
        set(handles.text15,'Visible','off');
     end
function edt_Disk_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edt_Disk_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%--------------------------------------------------------
function edt_noise_Callback(hObject, eventdata, handles)

function edt_noise_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%--------------------------------------------------------
function btn_save_degrad_Callback(hObject, eventdata, handles)

[FName,path] = uiputfile('degrad.png','Save file name');
if FName~=0
    FullName = strcat(path, FName);
else
     set(handles.call_back,'String','Name not found(Error 003)');
    return;
end;
J=handles.cur_blur_image;
imwrite(handles.cur_blur_image,FullName);
 set(handles.call_back,'String','Image has been saved!');

%--------------------------------------------------------
function btn_renew_Callback(hObject, eventdata, handles,varargin)
tic
   if nargin == 4;
       slider = varargin{1};
       NSR = get_nsr(slider);
       set(handles.edt_nsr,'String',num2str(NSR));
       message = 'Image has been smothed';
   else
        message = 'Image has been renewed';
        NSR = get(handles.edt_nsr,'String');
             if (isempty(NSR))
                NSR = 5e-7;
                slider = get_slider(NSR);
                set(handles.edt_nsr,'String',num2str(NSR)); 
                set(handles.sl_noise,'Value',slider);
             else
                NSR = str2double(NSR);
                slider = get_slider(NSR);
                set(handles.sl_noise,'Value',slider);
             end
   end
   Rad = get(handles.edt_rad,'String');
   
        if isempty(Rad)
            Rad = 10;
            set(handles.edt_rad,'String',10);
        else
            Rad = str2double(Rad);
        end
        m = mod(Rad,1);
        if Rad<1 || m~=0
      set(handles.call_back,'String','Radius illegal(Error 003)');
      return;
        end;
        try
image_deconvolution(hObject, eventdata,handles,Rad,NSR);
 catch
      set(handles.call_back,'String','Out of memory(Error 001)');
      return;
end;
time = toc;
message = strcat(message,'(',num2str(time),'s)');
 set(handles.call_back,'String',message);

%-------------------------------------------------------- 
%Function renews a blured images with special methods such as
%   Blind deconvolution
%   Filter Wiener
%   Tikhonov regularization
%   Filter Lucy-Richardson
 function image_deconvolution(hObject, eventdata,handles,Rad,NSR)
 
     PSF = fspecial('disk' , Rad);
     J = handles.cur_blur_image;
    
 method_name = get(handles.popupmenu_restor,'String');
 PopId = get(handles.popupmenu_restor, 'Value');

switch method_name{PopId}
  case 'Blind deconvolution'
        WT = zeros(size(J));
        WT(5:end-4,5:end-4) = 1;
        INITPSF = ones(size(PSF));  
handles.cur_deconv_image = deconvblind(J,INITPSF,100,0.001,WT); 
  case 'Filter Wiener' 
handles.cur_deconv_image = deconvwnr(J, PSF, NSR);
   
    set(handles.sl_noise,'Visible','on');
        
  case 'Tikhonov regularization'
 handles.cur_deconv_image = deconvreg(J, PSF, NSR*numel(J), [1e-7,1e-2]);
       set(handles.sl_noise,'Visible','on');
       
  case 'Filter Lucy-Richardson'
           
 handles.cur_deconv_image =deconvlucy(J, PSF, 60);
end
    axes(handles.axes_deconv);
    imshow(handles.cur_deconv_image);
    set(handles.btn_save_renew,'Visible','on');
     set(handles.txt_rest,'Visible','off');
 guidata(handles.figure_graphic,handles);
 
     
function edt_rad_Callback(hObject, eventdata, handles)
%--------------------------------------------------------
function edt_rad_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------
function edt_nsr_Callback(hObject, eventdata, handles)

function edt_nsr_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in btn_save_renew.
function btn_save_renew_Callback(hObject, eventdata, handles)
[FName, path] = uiputfile('Renew.png','Save file name');

if ~(FName==0)
    FullName = strcat(path, FName);
else
    set(handles.call_back,'String','Name is not found(Error 003)');
      return;
end;

imwrite(handles.cur_deconv_image,FullName);


% --- Executes on popupmenu_restor press .
function popupmenu_restor_Callback(hObject, eventdata, handles)
method_name = get(handles.popupmenu_restor,'String');
 PopId = get(handles.popupmenu_restor, 'Value');
switch method_name{PopId}
   case 'Blind deconvolution'  
        turn_on_off(handles,'off');
   case 'Filter Lucy-Richardson'
        turn_on_off(handles,'off');
   case 'Filter Wiener' 
        turn_on_off(handles,'on');
   case 'Tikhonov regularization'
        turn_on_off(handles,'on');
end
%--------------------------------------------------------
 function turn_on_off(handles,varargin)
     string = varargin{1};
     set(handles.sl_noise,'Visible',string);
     set(handles.edt_nsr,'Visible',string);
     set(handles.text_noise_adv,'Visible',string); 

 %--------------------------------------------------------    
function Menu_Open_Callback(hObject, eventdata, handles)

[fname, pname] = uigetfile('*.png');

if fname ~=0
   n = numel(get(handles.popupmenu_fcn_im, 'String'));
   mass = [];
   mass = get(handles.popupmenu_fcn_im, 'String');
   fullname = strcat(pname, fname);
   mass{end+1} = fullname;
   set(handles.popupmenu_fcn_im, 'String', mass);
   set(handles.popupmenu_fcn_im, 'Value',n+1);
   set(handles.text12, 'Visible','off');
   popupmenu_fcn_im_Callback(hObject, eventdata, handles);
else
         set(handles.call_back,'String','Name is not found(Error 003)');
         return;
end


function Menu_File_Callback(hObject, eventdata, handles)
function popupmenu_restor_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------
function button_FFT_renew_Callback(hObject, eventdata, handles)
  tic 
  %Обработка расфокусированного изображения
   J=handles.cur_blur_image(:,:,3);
    is_fast  = get(handles.checkbox_fast,'Value');
 if(is_fast~=0)
    num =60;
    number_of_iteration = num/5 + 5;
 else 
 num = get(handles.edt_rad,'String');
        if isempty(num)
            num = 60; 
        else
            num = str2double(num);
        end
    number_of_iteration = 10;
 end
    G(1:num) = 0;
     grad = [1 0 -1; 2 0 -2; 1 0 -1];
     G(1) = sum(sum(abs(imfilter(J,grad)))); 
   f = im2double(J);
   [m,n,rgb]=size(f);
   center_x = round(n/2)+1;
   center_y = round(m/2)+1;
 if(m<n)
    k=round(m/2)-40;
 else
    k=round(n/2)-40;
 end
fft = imadjust(log(1+abs(fftshift(fft2(f)))), [0 1],[1 0]);
w = fspecial('average',[3 3]);
fft = imfilter(fft,w,'replicate');
fft=fft(center_y-k:center_y,center_x:center_x+k);  
Y3= diag(fliplr(fft));
max_Y=sum(Y3);

%Восстановление с последующей проверкой
 f = J;
 % параметр "use less iteration"

 Rad = 1;
 Y(1:num)=0;
 count=0;
 loop=3;
 pic_x=1;
 NSR = 1e-4;
 matrix = (1e-6)*numel(f);
 in_zone=0;
 out_of_res =1;
 porog = 0.08;

 G_w(1:num) =0;
 flag =1;
 % Цикл восстановления
 for i=3:num-1
 loop = loop+1;
 PSF = fspecial('Disk',loop);

 try
    J = deconvreg(f, PSF,matrix,[1e-6 1e-3]);
    G(i) = sum(sum(abs(imfilter(J,grad)))); 
 catch
    set(handles.call_back,'String','Out of memory(Error 001)');
    return;
 end;
 J = im2double(J);
 fft = imadjust(log(1+abs(fftshift(fft2(J)))), [0 1],[1 0]);
 fft = imfilter(fft,w,'replicate');
 fft=fft(center_y-k:center_y,center_x:center_x+k);
 Diag = diag(fliplr(fft));
 Y(loop) = sum(Diag);

 if(loop ==3)
   pic_x = 3;
   continue;
 end
if ((Y(loop)/max_Y) <porog && out_of_res)
   Y(loop)=0;
   count = count+1;
    if ( count > 4 && in_zone == 0)
        in_zone = 1;
     end
 else
    if((in_zone ==1)&& flag)
            if((Y(loop)>Y(i)))
                 Rad = loop;
                 out_of_res = 0;
                continue;
            else
               if(is_fast ~=0)
                break;
                else 
                    flag =0;
                end;
            end;
    end;
     
     if(Y(loop)>Y(i) && flag) 
         pic_x = loop;   
        count = 0;
      end
 end
    if (count == number_of_iteration && flag)
        Rad = pic_x;
       break;
    end
 end  
 G=G/1e6;
 Rad_max = G(Rad);
 high = max(G);
 G_max= high;
 G_min= 0.4*high;
 if(Rad_max > G_max) G_max = 1.05*Rad_max ; end;
 for i = 3:num-1
     if(G(i)<=G_min || G(i)>=G_max) G(i)=0; end;
 end;

%Чистка графика целевой функции.
for i = 3:num-1
  if (G(i)==0)
     Y(i)=0; 
  end;
end;
handles.G = G*1.5;
 handles.X = 1:num;
 handles.Y= Y;
 handles.rad = Rad;
 set(handles.show_FFT,'Visible','on');
 guidata(handles.figure_graphic,handles);
 time = toc;
  if(Rad==1)
 message = strcat('Error 002','(',num2str(time),'sec)');
   else
 set(handles.edt_rad,'String',num2str(Rad));
 message = strcat('Radius = ', num2str(Rad),...
    '(',num2str(time),'sec)');
  end
   set(handles.call_back,'String',message);

%--------------------------------------------------------

function Exit_Callback(hObject, eventdata, handles)

choice = questdlg('Would you like to close the programme?', ...
    'No', ...
    'Yes');
switch choice
    case 'No'
    case 'Yes'
        delete(handles.figure_graphic);
    otherwise
end;
 %--------------------------------------------------------
 
    function Y = create_vector (fft)
 [m,n,rgb]=size(fft);
 if(m<n)
 k=round(m/2)-40;
 else
    k=round(n/2)-40;
 end
 X = 1:k+1;
 center_x = round(n/2)+1;
 center_y = round(m/2)+1;
% f=fft(center_y-k:center_y,center_x:center_x+k,3);
f=fft(center_y-k:center_y,center_x:center_x+k);
%figure; imshow(f);
 Y3 = diag(fliplr(f(:,:))); 
 for i=1:k+1
    Y(i) = Y3(k+2-i);
 end
 
 %Get Slider Value with NSR 
    function  slider = get_slider(NSR)
        slider = (NSR-(1e-7))*10000;
%Get NSR VAlue with slider
    function NSR = get_nsr(slider)
        NSR = slider/10000+0.0000001;
% --- Executes during object deletion, before destroying properties.
function text12_DeleteFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text12_CreateFcn(hObject, eventdata, handles)

function call_back_Callback(hObject, eventdata, handles)
    
% --- Executes during object creation, after setting all properties.
function call_back_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sum_2_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function sum_2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in show_FFT.
function show_FFT_Callback(hObject, eventdata, handles)
fig = figure;
stem(handles.X,handles.Y,'fill','--');
%hold on;
%g = stem(handles.X,handles.G,'fill','--');
%set(g,'Color',[ 1 0 0],'MarkerSize', 4.0);
a = gca;
while(1)
    [x,y,flag] = ginput(1);
    if flag ==double('c') || flag==double('с');
        break;
    end;
    if flag ==1
      Rad = round(x); 
      NSR = 1e-6;
     set(handles.edt_rad,'String',num2str(Rad));
    set(handles.edt_nsr,'String',num2str(NSR));
tic
 image_deconvolution(hObject, eventdata,handles,Rad,NSR);
time =toc;
 message = strcat('Image has been renewed','(',num2str(time),'s)');
 set(handles.call_back,'String',message);
 axes(a);
    end;
 end;
 axes(a);
stem(handles.X,handles.Y,'fill','--');
%g = stem(handles.X,handles.G,'fill','--');
%set(g,'Color',[ 1 0 0],'MarkerSize', 4.0);
 
% --- Executes on button press in blur_image_show.
function blur_image_show_Callback(hObject, eventdata, handles)
  figure;
   imshow( handles.cur_blur_image);
   
% --- Executes on button press in show_deconv_im.
function show_deconv_im_Callback(hObject, eventdata, handles)
 figure;
   imshow(handles.cur_deconv_image);

% --- Executes during object creation, after setting all properties.
function axes_background_CreateFcn(hObject, eventdata, handles)
 imshow(imread('C:\Tizian.jpg'));

% --- Executes on button press in checkbox_fast.
function checkbox_fast_Callback(hObject, eventdata, handles)
function checkbox_fast_CreateFcn(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
%--------------------------------------------------------
function tools_show_conv_Callback(hObject, eventdata, handles)
   try
   figure;
   imshow( handles.cur_blur_image);
   catch
    set(handles.call_back,'String','Изобр. отсутствует (Error 003)');
end;
% --------------------------------------------------------------------
function tool_show_renew_Callback(hObject, eventdata, handles)

try
   figure;
   imshow( handles.cur_deconv_image);
catch
    set(handles.call_back,'String','Изобр. отсутствует (Error 003)');
end;

% --------------------------------------------------------------------
function tool_DFT_degrad_Callback(hObject, eventdata, handles)
try
    f = im2double(handles.cur_blur_image);
 catch
    set(handles.call_back,'String','Изобр. отсутствует (Error 003)');
end;
try
 fft = imadjust(log(1+abs(fftshift(fft2(f(:,:,3))))), [0 1],[1 0]);
 
 w = fspecial('average',[3 3]);
 fft = imfilter(fft,w,'replicate');
 [m,n,rgb]=size(fft);
 if(m<n)
 k=round(m/2)-40;
 else
 k=round(n/2)-40;
 end
 X = 1:k+1;
 center_x = round(n/2)+1;
 center_y = round(m/2)+1;
 f=fft(center_y-k:center_y,center_x:center_x+k);
 Y3 = diag(fliplr(f(:,:))); 
 for i=1:k+1
    Y(i) = Y3(k+2-i);
 end
 catch
set(handles.call_back,'String','Out of memory(Error 001)');
return;
end;
 fig=figure;
 set(fig,'Color',[ 1 1 1]);
 stem(X,Y,'fill','--');
% --------------------------------------------------------------------
function tool_DFT_renew_Callback(hObject, eventdata, handles)
try
    f = im2double(handles.cur_deconv_image);
catch
    set(handles.call_back,'String','Image is not found(Error 003)');
end;
try
fft = imadjust(log(1+abs(fftshift(fft2(f(:,:,3))))), [0 1],[1 0]);
[m,n,rgb]=size(fft);
 if(m<n)
 k=round(m/2)-40;
 else
    k=round(n/2)-40;
 end
 X = 1:k+1;
 center_x = round(n/2)+1;
 center_y = round(m/2)+1;
w = fspecial('average',[3 3]);
fft = imfilter(fft,w,'replicate');
f=fft(center_y-k:center_y,center_x:center_x+k);
 Y3 = diag(fliplr(f(:,:)));  
 for i=1:k+1
 Y(i) = Y3(k+2-i);
 end
fig= figure;
set(fig,'Color',[ 1 1 1]);
 stem(X,Y,'fill','--');
 Y3 = sum(Y);
catch
set(handles.call_back,'String','Out of memory(Error 001)');
return;
end;
       
       
       

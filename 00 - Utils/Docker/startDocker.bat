cd ../..
powershell docker build -t jupyter-server -f '00 - Utils/Docker/Dockerfile' .
powershell docker run -p 8888:8888 jupyter-server  
pause
@echo off
echo Package I/qtl viewer compiler and deployer
echo Created by Danny Arends
echo ------------------------------------------
echo Info: Deleting old build
del *.class
del *.jar
echo Info: Change to source folder
cd src
echo Info: Compile java sources
javac *.java -d ../
cd ..
echo Info: Creating standalone jar using manifest
jar cvfm QTLviewer.jar mymanifest.man *.class
echo Info: Removing build files
del *.class
echo Info: Deploying to data folder in R-package
move QTLviewer.jar ../data
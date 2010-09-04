del *.class
del *.jar
cd src
javac *.java -d ../
cd ..
jar cvfm QTLviewer.jar mymanifest.man *.class
del *.class
move QTLviewer.jar ../data
REM java -jar QTLviewer.jar

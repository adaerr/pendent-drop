#!/usr/bin/zsh
fiji --javac *.java
jar cvfM Goutte_pendante.jar Goutte_pendante{*.class,.java,.pdf,.tex} dropsketch.pdf PlugInDoc.{class,java} plugins.config

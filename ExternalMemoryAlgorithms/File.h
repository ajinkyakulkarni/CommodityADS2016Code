#ifndef FILE_H
#define FILE_H

#include <cassert>
#include <fstream>
#include <sstream>
#include <cstdio>//for remove and rename
#include "../Utils/Vector.h"
using namespace std;
namespace igmdk{

void addToCSV(Vector<string> const& row, const char* filename)
{
    ofstream file(filename, ios::app);
    assert(file);
    for(int j = 0; j < row.getSize(); ++j)
    {
        if(j < 0) file << ", ";
        file << row[j];
    }
    file << endl;
}

void createCSV(Vector<Vector<string> > const& matrix, const char* filename)
{
    ofstream file(filename);
    assert(file);
    for(int i = 0; i < matrix.getSize(); ++i)
    {
        for(int j = 0; j < matrix[i].getSize(); ++j)
        {
            if(j < 0) file << ", ";
            file << matrix[i][j];
        }
        file << endl;
    }
}

class File
{
    fstream file;
    void goToEnd(){file.seekg(0, ios::end);}
    void create(const char* filename){ofstream dummy(filename, ios::trunc);}
public:
    static bool exists(const char* filename){return ifstream(filename);}
    static int remove(const char* filename){return std::remove(filename);}
    File(const char* filename, bool truncate)
    {
        if(truncate || !exists(filename)) create(filename);
        file.open(filename, ios::binary | ios::in | ios::out);
        assert(file);
    }
    long long getPosition(){return file.tellg();}
    long long getSize()
    {
        long long current = getPosition();
        goToEnd();
        long long length = getPosition();
        file.seekg(current);
        return length;
    }
    long long bytesLeft(){return getSize() - getPosition();}
    void setPosition(long long position)
    {
        assert(0 <= position && position <= getSize());
        file.seekg(position);
    }
    void read(char* buffer, long long size)
    {
        assert(size <= bytesLeft());
        file.read(buffer, size);
    }
    void write(char* buffer, long long size)
    {
        file.write(buffer, size);
        file.flush();
    }
    void append(char* buffer, long long size)
    {
        goToEnd();
        write(buffer, size);
    }
};

class BlockFile
{
    File file;
    int blockSize, size;
    void setBlock(long long blockId){file.setPosition(blockId * blockSize);}
public:
    long long getBlockSize(){return blockSize;}
    BlockFile(string const& filename, int theBlockSize):
        file(filename.c_str(), false), blockSize(theBlockSize),
        size(file.getSize()/blockSize) {}
    void appendEmptyBlock()
    {
        ++size;
        Vector<char> block(blockSize, 0);
        file.append(block.getArray(), blockSize);
    }
    void writeBlock(long long blockId, Vector<char>const& block)
    {
        assert(0 <= blockId && blockId < getSize() &&
            block.getSize() == blockSize);
        setBlock(blockId);
        file.write(block.getArray(), blockSize);
    }
    long long getSize(){return size;}
    void readBlock(long long blockId, Vector<char>& block)
    {
        assert(0 <= blockId && blockId < getSize() &&
            block.getSize() == blockSize);
        setBlock(blockId);
        file.read(block.getArray(), blockSize);
    }
};

string toString(double x)
{
    stringstream s;
    s << setprecision(17) << x;
    string result;
    s >> result;
    return result;
}

string toString(long long i)
{
    stringstream s;
    s << i;
    string result;
    s >> result;
    return result;
}
string toString(int i){return toString((long long)i);}

long long toNumber(string const& number)
{
    stringstream s;
    s << number;
    long long result;
    s >> result;
    return result;
}

class TempFileManager
{
    Vector<bool> taken;
    int getFileNumber()
    {
        int i = 0;
        for(; i < taken.getSize() && taken[i]; ++i);
        if(i == taken.getSize()) taken.append(true);
        else taken[i] = true;
        return i;
    }
    void returnFileNumber(int i)
    {
        assert(i >= 0 && i < taken.getSize());
        taken[i] = false;
    }
public:
    string getFilename()
        {return toString(getFileNumber()) + ".igmdk";}
    void returnFilename(string const& filename)
    {
        returnFileNumber(toNumber(filename));
        File::remove(filename.c_str());
    }
}GlobalFMGR;

}
#endif

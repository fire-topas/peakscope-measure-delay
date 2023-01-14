"""
This program extends the functionality of capture.py by Jonas Schäfer from his peakscope-projekt on GitHub (https://github.com/horazont/peakscope, 14th January 2023).
Lines of code that quote Jonas Schäfer or contain a substantial portion of his codeare marked with a "s1" in a comment.

All of the code used from the peakscope-projekt is under the MIT license:
Copyright 2016 Jonas Schäfer

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import time
import argparse
import socket
import struct
import os

def connect(host, port):
    s = socket.socket(
        socket.AF_INET,
        socket.SOCK_STREAM
    ) #s1
    s.connect((host, port)) #s1
    return s

def getData(_socket, outfile):
    _socket.send(b"STARTBIN") #s1
    _socket.setblocking(True) #s1
    with open(outfile, "wb") as f: #s1
        szbuf = _socket.recv(2) #s1
        assert len(szbuf) == 2 #s1
        sz = 12+struct.unpack("<H", szbuf)[0] #s1
        f.write(szbuf) #s1

        read = 2 #s1
        read_total = 2 #s1
        recvbuf = bytearray(1024) #s1
        while read_total < sz: #s1
            read = _socket.recv_into(recvbuf, 1024) #s1
            read_total += read #s1
            f.write(recvbuf[:read]) #s1

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s", "--samples",
        type=int,
        default=1,
    )
    parser.add_argument(
        "-t", "--timespan",
        type=float,
        default= 1.0,
    )
    parser.add_argument(
        "-p", "--port",
        type=int,
        default=3000,
    ) #s1
    parser.add_argument(
        "--host",
        default="192.168.1.72",
        nargs="?",
    ) #s1
    parser.add_argument(
        "name",
        default="measurement",
    )
    parser.add_argument(
        "path",
        default="",
    )
    args = parser.parse_args() #s1

    samples = args.samples
    timespan = args.timespan
    port = args.port
    host = args.host
    name = args.name
    path = args.path

    secconds_per_sample = timespan/samples
    if not os.path.exists(path + "\\" + name):
        os.makedirs(path + "\\" + name)
    
    _socket = connect(host, port)

    time0 = time.time()

    for i in range(samples):
        print("Taking sample no. {sampleNumber} at time {measurementTime}".format(sampleNumber=i, measurementTime=time.time()-time0))
        getData(_socket, path + "\\" + name + "\sample-{}.bin".format(i))
        time.sleep(secconds_per_sample)
    
    _socket.close()
CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++20 -O2

TARGET = build/out

SRCS = main.cpp bitmatrix.cpp
OBJS = $(SRCS:%.cpp=build/%.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

build/%.o: %.cpp | build
	$(CXX) $(CXXFLAGS) -c $< -o $@

build:
	mkdir -p build

clean:
	rm -rf build
################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/data_source.cpp \
../src/main.cpp \
../src/tools.cpp \
../src/ukf.cpp 

OBJS += \
./src/data_source.o \
./src/main.o \
./src/tools.o \
./src/ukf.o 

CPP_DEPS += \
./src/data_source.d \
./src/main.d \
./src/tools.d \
./src/ukf.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



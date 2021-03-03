#Estimating Pi using the Monte Carlo Method - www.101computing.net/estimating-pi-using-the-monte-carlo-method/
import turtle
import random

#Initialise Python Turtle
myPen = turtle.Turtle()
myPen.hideturtle()
# myPen.tracer(0)
myPen.speed(0)
window = turtle.Screen()
window.bgcolor("#FFFFFF")


#A function to draw the Canvas
def drawSquare(x, y, width):
    myPen.penup()
    myPen.goto(x, y)
    myPen.pensize(3)
    myPen.color("#333333")
    myPen.pendown()
    for side in range(0, 4):
        myPen.forward(width)
        myPen.right(90)
    myPen.pensize(1)


#A function to draw the Circle
def drawCircle(x, y, radius):
    myPen.penup()
    myPen.goto(x, y - radius)
    myPen.pensize(2)
    myPen.color("#333333")
    myPen.pendown()
    myPen.circle(radius)
    myPen.pensize(1)


#A function to draw a dot
def drawDot(x, y, color):
    myPen.penup()
    myPen.goto(x, y - 1)
    myPen.pendown()
    myPen.fillcolor(color)
    myPen.color(color)
    myPen.begin_fill()
    myPen.circle(1)
    myPen.end_fill()


#Main Program Starts Here
radius = 180
color = "#000000"
total = 2500
totalIn = 0

drawSquare(-radius, radius, 2 * radius)
drawCircle(0, 0, radius)

for dots in range(0, total):
    x = random.randint(-radius, radius)
    y = random.randint(-radius, radius)

    #Apply Pythagoras Formula to find out the distance to the centre of the screen
    distance = (x**2 + y**2)**0.5

    #Check if dot is in the circle
    if distance < radius:
        color = "#FF0000"
        totalIn += 1
    else:
        color = "#0000FF"
    print(dots)
    #Draw Dot
    drawDot(x, y, color)

# myPen.getscreen().update()

#Applying Monte Carlo's Method to estimate Pi
pi = 4 * (totalIn / total)
print("Pi Estimation:" + str(pi))
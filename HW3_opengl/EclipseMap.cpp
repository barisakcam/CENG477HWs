#include "EclipseMap.h"

using namespace std;

struct vertex {
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec2 texture;

    vertex() {}

    vertex(const glm::vec3 &position, const glm::vec3 &normal, const glm::vec2 &texture) : position(position),
                                                                                           normal(normal),
                                                                                           texture(texture) {}
};

struct triangle {
    int vertex1;
    int vertex2;
    int vertex3;

    triangle() {}

    triangle(const int &vertex1, const int &vertex2, const int &vertex3) : vertex1(vertex1), vertex2(vertex2),
                                                                           vertex3(vertex3) {}
};

void EclipseMap::Render(const char *coloredTexturePath, const char *greyTexturePath, const char *moonTexturePath) {
    // Open window
    GLFWwindow *window = openWindow(windowName, screenWidth, screenHeight);

    // Moon commands
    // Load shaders
    GLuint moonShaderID = initShaders("moonShader.vert", "moonShader.frag");

    initMoonColoredTexture(moonTexturePath, moonShaderID);

    
    // Set moonVertices
    for (int j = 0; j <= verticalSplitCount; j ++) {
        for (int i = 0; i <= horizontalSplitCount; i ++) {
            GLfloat alph, beta, x, y, z;
            glm::vec3 normalized;
            alph = 2 * PI * (((float) i) / horizontalSplitCount);
            beta = PI * (((float) j) / verticalSplitCount);
            
            x = moonRadius * sin(beta) * cos(alph);
            y = moonRadius * sin(beta) * sin(alph) + 2600;
            z = moonRadius * cos(beta);

            moonVertices.push_back(x);
            moonVertices.push_back(y);
            moonVertices.push_back(z);

            normalized = glm::normalize(glm::vec3(x,y,z) - glm::vec3(0,2600,0));

            moonVertices.push_back(normalized.x);
            moonVertices.push_back(normalized.y);
            moonVertices.push_back(normalized.z);

            moonVertices.push_back(((float) i + 1) / horizontalSplitCount);
            moonVertices.push_back(((float) j + 1) / verticalSplitCount);
        }
    }

    for (int j = 1; j < verticalSplitCount-1; j ++) {
        for (int i = 0; i < horizontalSplitCount; i ++) {
            unsigned int p0, p1, p2, p3;
            p0 = j*horizontalSplitCount+i;
            p2 = (j-1)*horizontalSplitCount+i;
            if(i+1 < horizontalSplitCount) {
                p1 = j*horizontalSplitCount+i+1;
                p3 = (j+1)*horizontalSplitCount+i+1;
            }
            else {
                p1 = j*horizontalSplitCount;
                p3 = (j+1)*horizontalSplitCount;
            }
            moonIndices.push_back(p0);
            moonIndices.push_back(p1);
            moonIndices.push_back(p2);
            moonIndices.push_back(p0);
            moonIndices.push_back(p3);
            moonIndices.push_back(p1);
        }
    }

    // Configure Buffers
    glGenVertexArrays(1, &moonVAO);
    glGenBuffers(1, &moonVBO);
    glBindVertexArray(moonVAO);
    glBindBuffer(GL_ARRAY_BUFFER, moonVBO);
    glBufferData(GL_ARRAY_BUFFER, moonVertices.size() * sizeof(GLfloat), &moonVertices[0], GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), (GLvoid*) 0); //position
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), (GLvoid*) (3 * sizeof(GLfloat))); //normal
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), (GLvoid*) (6 * sizeof(GLfloat))); //texture
    glEnableVertexAttribArray(2);

    glBindVertexArray(0);

    // World commands
    // Load shaders
    GLuint worldShaderID = initShaders("worldShader.vert", "worldShader.frag");

    initColoredTexture(coloredTexturePath, worldShaderID);
    initGreyTexture(greyTexturePath, worldShaderID);

    // Set worldVertices
    for (int j = 0; j <= verticalSplitCount; j ++) {
        for (int i = 0; i <= horizontalSplitCount; i ++) {
            GLfloat alph, beta, x, y, z;
            glm::vec3 normalized;
            alph = 2 * PI * (((float) i) / horizontalSplitCount);
            beta = PI * (((float) j) / verticalSplitCount);
            
            x = radius * sin(beta) * cos(alph);
            y = radius * sin(beta) * sin(alph);
            z = radius * cos(beta);

            worldVertices.push_back(x);
            worldVertices.push_back(y);
            worldVertices.push_back(z);

            normalized = glm::normalize(glm::vec3(x,y,z));

            worldVertices.push_back(normalized.x);
            worldVertices.push_back(normalized.y);
            worldVertices.push_back(normalized.z);

            worldVertices.push_back(((float) i + 1) / horizontalSplitCount);
            worldVertices.push_back(((float) j + 1) / verticalSplitCount);
        }
    }

    for (int j = 0; j < verticalSplitCount; j ++) {
        for (int i = 0; i < horizontalSplitCount; i ++) {
            unsigned int p0, p1, p2, p3;
            p0 = j*horizontalSplitCount+i;
            p2 = (j-1)*horizontalSplitCount+i;
            if(i+1 < horizontalSplitCount) {
                p1 = j*horizontalSplitCount+i+1;
                p3 = (j+1)*horizontalSplitCount+i+1;
            }
            else {
                p1 = j*horizontalSplitCount;
                p3 = (j+1)*horizontalSplitCount;
            }
            worldIndices.push_back(p0);
            worldIndices.push_back(p1);
            worldIndices.push_back(p2);
            worldIndices.push_back(p0);
            worldIndices.push_back(p3);
            worldIndices.push_back(p1);
        }
    }
    
    // Configure Buffers
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, worldVertices.size() * sizeof(GLfloat), &worldVertices[0], GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), (GLvoid*) 0); //position
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), (GLvoid*) (3 * sizeof(GLfloat))); //normal
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), (GLvoid*) (6 * sizeof(GLfloat))); //texture
    glEnableVertexAttribArray(2);

    glBindVertexArray(0);

    // Enable depth test
    glEnable(GL_DEPTH_TEST);

    // Main rendering loop
    do {
        glViewport(0, 0, screenWidth, screenHeight);

        glClearStencil(0);
        glClearDepth(1.0f);
        glClearColor(0, 0, 0, 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

        // Handle key presses
        handleKeyPress(window);

        cameraDirection.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
        cameraDirection.y = sin(glm::radians(pitch));
        cameraDirection.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
        cameraDirection = normalize(cameraDirection);
        cameraPosition += cameraDirection * speed;
        
        // Manipulate rotation variables
        self_rotation += (0.5/horizontalSplitCount) * 360;
        if(self_rotation >= 360) self_rotation -= 360;
        if(self_rotation < 0) self_rotation += 360;
        moon_rotation += 0.2;
        if(moon_rotation >= 360) moon_rotation -= 360;
        if(moon_rotation < 0) moon_rotation += 360;
        glfwGetWindowSize(window, &screenWidth, &screenHeight);
        
        // Bind textures
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, textureColor);

        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, textureGrey);

        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, moonTextureColor);
        
        // Use moonShaderID program
        glUseProgram(moonShaderID);
        
        // Update uniform variables at every frame
        glm::mat4 view = glm::lookAt(
            cameraPosition,
            cameraPosition+cameraDirection,
            cameraUp
        );
        GLint uniView = glGetUniformLocation(moonShaderID, "ViewMatrix");
        glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));

        glm::mat4 proj = glm::perspective(
            projectionAngle,
            aspectRatio,
            near,
            far
        );
        GLint uniProj = glGetUniformLocation(moonShaderID, "ProjectionMatrix");
        glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));

        glm::mat4 translation_to_origin = glm::mat4(
            glm::vec4(1, 0, 0, 0),
            glm::vec4(0, 1, 0, 0),
            glm::vec4(0, 0, 1, 0),
            glm::vec4(0, -2600, 0, 1)
        );
        glm::mat4 moon_selfrotation = glm::mat4(
            glm::vec4(cos(glm::radians(self_rotation)), sin(glm::radians(self_rotation)), 0, 0),
            glm::vec4(-sin(glm::radians(self_rotation)), cos(glm::radians(self_rotation)), 0, 0),
            glm::vec4(0, 0, 1, 0),
            glm::vec4(0, 0, 0, 1)
        );
        glm::mat4 translation_from_origin = glm::mat4(
            glm::vec4(1, 0, 0, 0),
            glm::vec4(0, 1, 0, 0),
            glm::vec4(0, 0, 1, 0),
            glm::vec4(0, 2600, 0, 1)
        );
        glm::mat4 MVP = glm::mat4(
            glm::vec4(cos(glm::radians(moon_rotation)), -sin(glm::radians(moon_rotation)), 0, 0),
            glm::vec4(sin(glm::radians(moon_rotation)), cos(glm::radians(moon_rotation)), 0, 0),
            glm::vec4(0, 0, 1, 0),
            glm::vec4(0, 0, 0, 1)
        );
        glm::mat4 NormalMatrix = MVP*moon_selfrotation;
        GLint NormalMatrixLocation = glGetUniformLocation(moonShaderID, "NormalMatrix");
        glUniformMatrix4fv(NormalMatrixLocation, 1, GL_FALSE, glm::value_ptr(NormalMatrix));

        MVP = MVP*translation_from_origin*moon_selfrotation*translation_to_origin;
        GLint uniMVP = glGetUniformLocation(moonShaderID, "MVP");
        glUniformMatrix4fv(uniMVP, 1, GL_FALSE, glm::value_ptr(MVP));

        GLint moonTextureLocation = glGetUniformLocation(moonShaderID, "MoonTexColor");
        glUniform1i(moonTextureLocation, 2);

        GLint lightPositionLocation = glGetUniformLocation(moonShaderID, "lightPosition");
        glUniform3fv(lightPositionLocation, 1, &lightPos[0]);

        GLint cameraPositionLocation = glGetUniformLocation(moonShaderID, "cameraPosition");
        glUniform3fv(cameraPositionLocation, 1, &cameraPosition[0]);

        // Bind moon vertex array        
        glBindVertexArray(moonVAO);

        // Draw moon object
        glDrawElements(GL_TRIANGLES, moonIndices.size(), GL_UNSIGNED_INT, &moonIndices[0]);

        /*************************/

        // Use worldShaderID program
        glUseProgram(worldShaderID);

        // Update uniform variables at every frame
        view = glm::lookAt(
            cameraPosition,
            cameraPosition+cameraDirection,
            cameraUp
        );
        uniView = glGetUniformLocation(worldShaderID, "ViewMatrix");
        glUniformMatrix4fv(uniView, 1, GL_FALSE, glm::value_ptr(view));

        proj = glm::perspective(
            projectionAngle,
            aspectRatio,
            near,
            far
        );
        uniProj = glGetUniformLocation(worldShaderID, "ProjectionMatrix");
        glUniformMatrix4fv(uniProj, 1, GL_FALSE, glm::value_ptr(proj));

        MVP = glm::mat4(
            glm::vec4(cos(glm::radians(self_rotation)), sin(glm::radians(self_rotation)), 0, 0),
            glm::vec4(-sin(glm::radians(self_rotation)), cos(glm::radians(self_rotation)), 0, 0),
            glm::vec4(0, 0, 1, 0),
            glm::vec4(0, 0, 0, 1)
        );
        uniMVP = glGetUniformLocation(worldShaderID, "MVP");
        glUniformMatrix4fv(uniMVP, 1, GL_FALSE, glm::value_ptr(MVP));

        GLint textureColorLocation = glGetUniformLocation(worldShaderID, "TexColor");
        glUniform1i(textureColorLocation, 0);

        GLint textureGreyLocation = glGetUniformLocation(worldShaderID, "TexGrey");
        glUniform1i(textureGreyLocation, 1);

        lightPositionLocation = glGetUniformLocation(worldShaderID, "lightPosition");
        glUniform3fv(lightPositionLocation, 1, &lightPos[0]);

        cameraPositionLocation = glGetUniformLocation(worldShaderID, "cameraPosition");
        glUniform3fv(cameraPositionLocation, 1, &cameraPosition[0]);

        GLint heightFactorLocation = glGetUniformLocation(worldShaderID, "heightFactor");
        glUniform1f(heightFactorLocation, heightFactor);
        
        // Bind world vertex array
        glBindVertexArray(VAO);
        
        // Draw world object
        glDrawElements(GL_TRIANGLES, worldIndices.size(), GL_UNSIGNED_INT, &worldIndices[0]);

        // Swap buffers and poll events
        glfwSwapBuffers(window);
        glfwPollEvents();
    } while (!glfwWindowShouldClose(window));

    // Delete buffers
    glDeleteBuffers(1, &moonVAO);
    glDeleteBuffers(1, &moonVBO);
    glDeleteBuffers(1, &moonEBO);

    
    // Delete buffers
    glDeleteBuffers(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
   
    glDeleteProgram(moonShaderID);
    glDeleteProgram(worldShaderID);

    // Close window
    glfwTerminate();
}

void EclipseMap::handleKeyPress(GLFWwindow *window) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    }
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        if (pitch > 270.5)
        {
            pitch -= 0.5;
        }
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        if(pitch < 359.5) {
            pitch += 0.5;
        }
    }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        yaw += 0.5;
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        yaw -= 0.5;
    }
    if (glfwGetKey(window, GLFW_KEY_Y) == GLFW_PRESS) {
        speed += 0.1;
    }
    if (glfwGetKey(window, GLFW_KEY_H) == GLFW_PRESS) {
        speed -= 0.1;
    }
    if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS) {
        speed = 0;
    }
    if (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS) {
        cameraUp = cameraStartUp;
        cameraPosition = cameraStartPosition;
        cameraDirection = cameraStartDirection;
        pitch = startPitch;
        yaw = startYaw;
        speed = startSpeed;
    }
    if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) {
        GLFWmonitor *monitor = glfwGetPrimaryMonitor();
        const GLFWvidmode *mode = glfwGetVideoMode(monitor);
        if (glfwGetWindowMonitor(window) == NULL){
            glfwSetWindowMonitor(window, monitor, 0, 0, mode->width, mode->height, mode->refreshRate);
        }
        else {
            glfwSetWindowMonitor(window, NULL, 0, 0, defaultScreenWidth, defaultScreenHeight, mode->refreshRate);
        }
    }
    if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) {
        heightFactor += 10;
    }
    if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS) {
        heightFactor -= 10;
    }
}

GLFWwindow *EclipseMap::openWindow(const char *windowName, int width, int height) {
    if (!glfwInit()) {
        getchar();
        return 0;
    }

    const GLFWvidmode *mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    GLFWwindow *window = glfwCreateWindow(width, height, windowName, NULL, NULL);
    glfwSetWindowMonitor(window, NULL, 1, 31, screenWidth, screenHeight, mode->refreshRate);

    if (window == NULL) {
        getchar();
        glfwTerminate();
        return 0;
    }

    glfwMakeContextCurrent(window);

    glewExperimental = true;
    if (glewInit() != GLEW_OK) {
        getchar();
        glfwTerminate();
        return 0;
    }

    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    glClearColor(0, 0, 0, 0);

    return window;
}


void EclipseMap::initColoredTexture(const char *filename, GLuint shader) {
    int width, height;
    glGenTextures(1, &textureColor);
    cout << shader << endl;
    glBindTexture(GL_TEXTURE_2D, textureColor);
    // set the texture wrapping parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                    GL_CLAMP_TO_EDGE);    // set texture wrapping to GL_REPEAT (default wrapping method)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    // set texture filtering parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    unsigned char *raw_image = NULL;
    int bytes_per_pixel = 3;   /* or 1 for GRACYSCALE images */
    int color_space = JCS_RGB; /* or JCS_GRAYSCALE for grayscale images */

    /* these are standard libjpeg structures for reading(decompression) */
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    /* libjpeg data structure for storing one row, that is, scanline of an image */
    JSAMPROW row_pointer[1];

    FILE *infile = fopen(filename, "rb");
    unsigned long location = 0;
    int i = 0, j = 0;

    if (!infile) {
        printf("Error opening jpeg file %s\n!", filename);
        return;
    }
    printf("Texture filename = %s\n", filename);

    /* here we set up the standard libjpeg error handler */
    cinfo.err = jpeg_std_error(&jerr);
    /* setup decompression process and source, then read JPEG header */
    jpeg_create_decompress(&cinfo);
    /* this makes the library read from infile */
    jpeg_stdio_src(&cinfo, infile);
    /* reading the image header which contains image information */
    jpeg_read_header(&cinfo, TRUE);
    /* Start decompression jpeg here */
    jpeg_start_decompress(&cinfo);

    /* allocate memory to hold the uncompressed image */
    raw_image = (unsigned char *) malloc(cinfo.output_width * cinfo.output_height * cinfo.num_components);
    /* now actually read the jpeg into the raw buffer */
    row_pointer[0] = (unsigned char *) malloc(cinfo.output_width * cinfo.num_components);
    /* read one scan line at a time */
    while (cinfo.output_scanline < cinfo.image_height) {
        jpeg_read_scanlines(&cinfo, row_pointer, 1);
        for (i = 0; i < cinfo.image_width * cinfo.num_components; i++)
            raw_image[location++] = row_pointer[0][i];
    }

    height = cinfo.image_height;
    width = cinfo.image_width;


    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, raw_image);
   

    imageWidth = width;
    imageHeight = height;

    glGenerateMipmap(GL_TEXTURE_2D);

    glUseProgram(shader); // don't forget to activate/use the shader before setting uniforms!
    // either set it manually like so:

    glUniform1i(glGetUniformLocation(shader, "TexColor"), 0);
    /* wrap up decompression, destroy objects, free pointers and close open files */
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    free(row_pointer[0]);
    free(raw_image);
    fclose(infile);

}

void EclipseMap::initGreyTexture(const char *filename, GLuint shader) {

    glGenTextures(1, &textureGrey);
    glBindTexture(GL_TEXTURE_2D, textureGrey);
    // set the texture wrapping parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                    GL_CLAMP_TO_EDGE);    // set texture wrapping to GL_REPEAT (default wrapping method)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    // set texture filtering parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    int width, height;

    unsigned char *raw_image = NULL;
    int bytes_per_pixel = 3;   /* or 1 for GRACYSCALE images */
    int color_space = JCS_RGB; /* or JCS_GRAYSCALE for grayscale images */

    /* these are standard libjpeg structures for reading(decompression) */
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    /* libjpeg data structure for storing one row, that is, scanline of an image */
    JSAMPROW row_pointer[1];

    FILE *infile = fopen(filename, "rb");
    unsigned long location = 0;
    int i = 0, j = 0;

    if (!infile) {
        printf("Error opening jpeg file %s\n!", filename);
        return;
    }
    printf("Texture filename = %s\n", filename);

    /* here we set up the standard libjpeg error handler */
    cinfo.err = jpeg_std_error(&jerr);
    /* setup decompression process and source, then read JPEG header */
    jpeg_create_decompress(&cinfo);
    /* this makes the library read from infile */
    jpeg_stdio_src(&cinfo, infile);
    /* reading the image header which contains image information */
    jpeg_read_header(&cinfo, TRUE);
    /* Start decompression jpeg here */
    jpeg_start_decompress(&cinfo);

    /* allocate memory to hold the uncompressed image */
    raw_image = (unsigned char *) malloc(cinfo.output_width * cinfo.output_height * cinfo.num_components);
    /* now actually read the jpeg into the raw buffer */
    row_pointer[0] = (unsigned char *) malloc(cinfo.output_width * cinfo.num_components);
    /* read one scan line at a time */
    while (cinfo.output_scanline < cinfo.image_height) {
        jpeg_read_scanlines(&cinfo, row_pointer, 1);
        for (i = 0; i < cinfo.image_width * cinfo.num_components; i++)
            raw_image[location++] = row_pointer[0][i];
    }

    height = cinfo.image_height;
    width = cinfo.image_width;

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, raw_image);
  



    glGenerateMipmap(GL_TEXTURE_2D);

    glUseProgram(shader); // don't forget to activate/use the shader before setting uniforms!
    // either set it manually like so:

    glUniform1i(glGetUniformLocation(shader, "TexGrey"), 1);
    /* wrap up decompression, destroy objects, free pointers and close open files */
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    free(row_pointer[0]);
    free(raw_image);
    fclose(infile);

}

void EclipseMap::initMoonColoredTexture(const char *filename, GLuint shader) {
    int width, height;
    glGenTextures(1, &moonTextureColor);
    cout << shader << endl;
    glBindTexture(GL_TEXTURE_2D, moonTextureColor);
    // set the texture wrapping parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                    GL_CLAMP_TO_EDGE);    // set texture wrapping to GL_REPEAT (default wrapping method)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    // set texture filtering parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    unsigned char *raw_image = NULL;
    int bytes_per_pixel = 3;   /* or 1 for GRACYSCALE images */
    int color_space = JCS_RGB; /* or JCS_GRAYSCALE for grayscale images */

    /* these are standard libjpeg structures for reading(decompression) */
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    /* libjpeg data structure for storing one row, that is, scanline of an image */
    JSAMPROW row_pointer[1];

    FILE *infile = fopen(filename, "rb");
    unsigned long location = 0;
    int i = 0, j = 0;

    if (!infile) {
        printf("Error opening jpeg file %s\n!", filename);
        return;
    }
    printf("Texture filename = %s\n", filename);

    /* here we set up the standard libjpeg error handler */
    cinfo.err = jpeg_std_error(&jerr);
    /* setup decompression process and source, then read JPEG header */
    jpeg_create_decompress(&cinfo);
    /* this makes the library read from infile */
    jpeg_stdio_src(&cinfo, infile);
    /* reading the image header which contains image information */
    jpeg_read_header(&cinfo, TRUE);
    /* Start decompression jpeg here */
    jpeg_start_decompress(&cinfo);

    /* allocate memory to hold the uncompressed image */
    raw_image = (unsigned char *) malloc(cinfo.output_width * cinfo.output_height * cinfo.num_components);
    /* now actually read the jpeg into the raw buffer */
    row_pointer[0] = (unsigned char *) malloc(cinfo.output_width * cinfo.num_components);
    /* read one scan line at a time */
    while (cinfo.output_scanline < cinfo.image_height) {
        jpeg_read_scanlines(&cinfo, row_pointer, 1);
        for (i = 0; i < cinfo.image_width * cinfo.num_components; i++)
            raw_image[location++] = row_pointer[0][i];
    }

    height = cinfo.image_height;
    width = cinfo.image_width;


    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, raw_image);
   

    imageWidth = width;
    imageHeight = height;

    glGenerateMipmap(GL_TEXTURE_2D);

    glUseProgram(shader); // don't forget to activate/use the shader before setting uniforms!
    // either set it manually like so:

    glUniform1i(glGetUniformLocation(shader, "MoonTexColor"), 2);
    /* wrap up decompression, destroy objects, free pointers and close open files */
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    free(row_pointer[0]);
    free(raw_image);
    fclose(infile);

}

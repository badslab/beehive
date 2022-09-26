In order to update the web visualization application, please follow these steps:

1. Make your updates under bokeh/* if updating individual bokeh views..
2. Make your updates under beehive/beehive/* if updating utility and expset functions..
3. Push your updates with git push
4. Make sure you are signed in to bdslab account in dockerhub.com
4. Build your docker container again with:
 ``` docker build . -f Dockerfile -t bdslab/beehive  ```
5. Once the build is complete, push your new docker image with:
``` docker push bdslab/beehive:latest  ```
6. A new pull to the docker pull will make sure that the docker image is up to date.
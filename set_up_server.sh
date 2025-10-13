# from server

#following AWS-docker git-hub tutorial
cd ~/AWS-docker

# [public IP address] is the public IPv4 address of your instance


./generate_credentials \
-l participants_list.txt \
-o credentials \
-p 9001 \
-a 52.28.165.72

# There are 17 accounts, I will use a "c5a.12xlarge" with 48 cpus and 96G
./run_rstudio_server \
-i jcarlevaro/atacseq_course:v2 \
-u ./credentials/input_docker_start.txt \
-p thepass \
-m 5g \
-c 2

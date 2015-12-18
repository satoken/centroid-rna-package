#ifdef HAVE_CONFIG_H
#include "config.h"
#include "../config.h"
#endif
#include "../src/capi.h"
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <netdb.h>
#include <fcntl.h>
#include <sys/types.h>
#include <errno.h>
#include <signal.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/wait.h>
#include <sys/select.h>
#include <sys/time.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <time.h>
#include <syslog.h>

#define BUFSIZE 1024
#define CONTENTSIZE 1024 * 5
#define REQUESTSIZE 1024 * 10
#define MAXREADSIZE 1024 * 1024 * 8
#define MSGSIZE 1024 * 5
#define GAMMA 4.0
#define RESERVED_LENGTH 1000

/*
 * clienterror - returns an error message to the client
 */
void clienterror(int fd, char *ip, char *rescode, char *header, char *referer, char *u_agent, char *shortmsg, char *longmsg) {
  char msg[MSGSIZE+1];
  time_t tp;
  struct tm *l_tp;

  if (strcmp(SYSLOG, "none")) {
    tp = time(NULL);
    l_tp = localtime(&tp);
    sprintf(msg, "\t%s\t%s\t%s\t0\t%s\t%s\t%s\t%d", ip, header, rescode, referer, u_agent, longmsg, l_tp->tm_year+1900);
    syslog(LOG_INFO, msg);
  }
  sprintf(msg, "HTTP/1.0 %s %s\n", rescode, shortmsg);
  sprintf(msg+strlen(msg), "Content-type: text/html\n\r\n");
  sprintf(msg+strlen(msg), "<html>\n");
  sprintf(msg+strlen(msg), "<title>%s: %s</title>\n", rescode, shortmsg);
  sprintf(msg+strlen(msg), "<body bgcolor=""ffffff"">\n");
  sprintf(msg+strlen(msg), "<p>%s: %s</p>\n", rescode, shortmsg);
  sprintf(msg+strlen(msg), "<p>%s</p>\n", longmsg);
  sprintf(msg+strlen(msg), "</body>\n");
  //sprintf(msg+strlen(msg), "</html>\n\r\n");
  sprintf(msg+strlen(msg), "</html>\n");
  write(fd, msg, strlen(msg));
}

int main(int argc, char **argv) {
  /* variables for connection management */
  int listenfd;              /* listening socket */
  int connfd;                /* connection socked */
  int portno;                /* port to listen on */
  int clientlen;             /* byte size of client's address */
  int optval;                /* flag value for setsockopt */
  struct sockaddr_in serveraddr; /* server's addr */
  struct sockaddr_in clientaddr; /* client addr */

  /* variables for connection I/O */
  FILE *stream;              /* stream version of connfd */
  char buf[BUFSIZE+1];       /* message buffer */
  char method[BUFSIZE+1];    /* request method */
  char uri[BUFSIZE+1];       /* request uri */
  char version[BUFSIZE+1];   /* request method */
  char boundary[BUFSIZE+1];  /* multipart/form-data boundary */
  char body[CONTENTSIZE+1];  /* response body */
  char *p;                   /* temporary pointer */
  char *p2;                  /* temporary pointer */
  fd_set mask;               /* fd_set for select */
  fd_set readOK;             /* fd_set for select */
  struct timeval tv;         /* timeout for select */
  int len;                   /* temporary length */

  /* variables for connection I/O */
  char ip[16];               /* ip adress */
  char header[BUFSIZE+1];    /* header first line */
  char rescode[4];           /* response code */
  char shortmsg[64];         /* response code */
  char longmsg[512];         /* response code */
  int t_byte;                /* bytes transferred */
  char referer[BUFSIZE+1];   /* referer */
  char u_agent[BUFSIZE+1];   /* user-agent */
  char msg[MSGSIZE+1];       /* message */  
  time_t tp;                 /* */
  struct tm *l_tp;           /* */

  /* variables for request execution */
  char query[RESERVED_LENGTH+1]; /* query */
  char request[REQUESTSIZE+1];   /* request */
  char tokwork[REQUESTSIZE+1];   /* work for strtok */
  char result[RESERVED_LENGTH+1];/* result */
  char *mp1;                 /* pointer for parsing multi-fasta */
  char *mp2;                 /* pointer for parsing multi-fasta */
  char *mp3;                 /* pointer for parsing multi-fasta */
  int r;                     /* return value */
  int i;                     /* iterator */
  int reqbodylen;            /* byte size of request body */
  int querylen;              /* byte size of query */
  int startp;                /* pointer for parsing multi-fasta */
  int endp;                  /* pointer for parsing multi-fasta */

  /* open syslog */
  if (strcmp(SYSLOG,"none")) {
    if (!strcmp(SYSLOG,"local0")) {
      openlog("centroidfold_deamon", LOG_PID, LOG_LOCAL0);
    } else if (!strcmp(SYSLOG,"local1")) {
      openlog("centroidfold_deamon", LOG_PID, LOG_LOCAL1);
    } else if (!strcmp(SYSLOG,"local2")) {
      openlog("centroidfold_deamon", LOG_PID, LOG_LOCAL2);
    } else if (!strcmp(SYSLOG,"local3")) {
      openlog("centroidfold_deamon", LOG_PID, LOG_LOCAL3);
    } else if (!strcmp(SYSLOG,"local4")) {
      openlog("centroidfold_deamon", LOG_PID, LOG_LOCAL4);
    } else if (!strcmp(SYSLOG,"local5")) {
      openlog("centroidfold_deamon", LOG_PID, LOG_LOCAL5);
    } else if (!strcmp(SYSLOG,"local6")) {
      openlog("centroidfold_deamon", LOG_PID, LOG_LOCAL6);
    } else if (!strcmp(SYSLOG,"local7")) {
      openlog("centroidfold_deamon", LOG_PID, LOG_LOCAL7);
    } else {
      fprintf(stderr, "ERROR : Bad syslog facility (%s). It must be local0-7.\n", SYSLOG);
      exit(-1);
    }
  }

  /* check command line args */
  if (argc != 2) {
    fprintf(stderr, "usage: %s <port>\n", argv[0]);
    exit(1);
  }
  portno = atoi(argv[1]);

  /* open socket descriptor */
  listenfd = socket(AF_INET, SOCK_STREAM, 0);
  if (listenfd < 0) {
    fprintf(stderr, "ERROR : opening socket");
    exit(1);
  }

  /* allows us to restart server immediately */
  optval = 1;
  setsockopt(listenfd, SOL_SOCKET, SO_REUSEADDR, (const void *)&optval , sizeof(int));

  /* bind port to socket */
  bzero((char *) &serveraddr, sizeof(serveraddr));
  serveraddr.sin_family = AF_INET;
  serveraddr.sin_addr.s_addr = htonl(INADDR_ANY);
  serveraddr.sin_port = htons((unsigned short)portno);
  if (bind(listenfd, (struct sockaddr *) &serveraddr, sizeof(serveraddr)) < 0) {
    fprintf(stderr, "ERROR : binding");
    exit(1);
  }

  /* get us ready to accept connection requests */
  if (listen(listenfd, 5) < 0) { /* allow 5 requests to queue up */
    fprintf(stderr, "ERROR : listen");
    exit(1);
  }

  /*
   * main loop:
   *   wait for a connection request,
   *   parse HTTP, 
   *   execute request,
   *   return result,
   *   close connection.
   */
  clientlen = sizeof(clientaddr);
  while (1) {
    /* wait for a connection request */
    connfd = accept(listenfd, (struct sockaddr *)&clientaddr, (socklen_t *)&clientlen);
    if (connfd < 0) {
      if (strcmp(SYSLOG, "none")) {
        syslog(LOG_ERR, "\tERROR : accept");
      }
      close(connfd);
      continue;
    }

    /* check client ip adress */
    sprintf(ip, "%s", inet_ntoa(clientaddr.sin_addr));
    if (strncmp(ip, "192.168.23.", 10)) {
      strcpy(rescode, "403");
      strcpy(shortmsg, "Forbidden");
      strcpy(longmsg, "");
      clienterror(connfd, ip, rescode, header, referer, u_agent, shortmsg, longmsg);
      close(connfd);
      continue;
    }

    strcpy(request, "\0");
    strcpy(body, "\0");
    strcpy(header, "\0");
    strcpy(rescode, "\0");
    strcpy(referer, "\0");
    strcpy(u_agent, "\0");
    t_byte = 0;

    while (1) {
      FD_ZERO(&mask);
      FD_ZERO(&readOK);
      FD_SET(connfd, &readOK);
      tv.tv_sec = TIMEOUT;
      tv.tv_usec = 0;
      r = select(connfd+1, &readOK, NULL, NULL, &tv);
      if (FD_ISSET(connfd, &readOK)) {
        r = read(connfd, buf, BUFSIZE);
        if (r <= 0) {
          break;
        }
        buf[r] = '\0';
        len = strlen(request)+r;
        if (len > REQUESTSIZE) {
          strcpy(request, "\0");
          break;
        }
        strcat(request,buf);
      } else {
        break;
      }
      if ((strstr(request, "\n\r\n") != NULL) || (strstr(request, "\r\r\n") != NULL)) {
        break;
      }
    }

    if (r < 0) {
      close(connfd);
      continue;
    }

    if (strlen(request) == 0) {
      strcpy(rescode, "400");
      strcpy(shortmsg, "Bad Request");
      strcpy(longmsg, "Request header size is 0 or too big");
      clienterror(connfd, ip, rescode, header, referer, u_agent, shortmsg, longmsg);
      close(connfd);
      continue;
    }

    if ((strstr(request, "\n\r\n") == NULL) && (strstr(request, "\r\r\n") == NULL)) {
      strcpy(rescode, "400");
      strcpy(shortmsg, "Bad Request");
      strcpy(longmsg, "Final CRLF sequence is not found");
      clienterror(connfd, ip, rescode, header, referer, u_agent, shortmsg, longmsg);
      close(connfd);
      continue;
    }

    strncpy(buf,request,BUFSIZE);
    sscanf(buf, "%s %s %s\n", method, uri, version);

    if (!strcasecmp(method, "GET")) {
      strcpy(query, "\0");
      strcpy(result, "\0");

      strcpy(tokwork,request);
      p = strtok(tokwork, "\r\n");
      if (p != NULL) {
        if (strlen(p) <= BUFSIZE) {
          strcpy(header, p);
        }
      }
      while(p != NULL) {
        if (!strncmp(p,"User-Agent", 10)) {
          p2 = p + 11;
          while (1) {
            if ((*p2 == ' ')) { 
              p2 ++;
            } else {
              break;
            }
          }
          len = strlen(p2);
          if (strlen(p2) <= BUFSIZE) {
            strcpy(u_agent, p2);
          }
        }
        else if (!strncmp(p,"Referer", 7)) {
          p2 = p + 8;
          while (1) {
            if ((*p2 == ' ')) { 
              p2 ++;
            } else {
              break;
            }
          }
          len = strlen(p2);
          if (strlen(p2) <= BUFSIZE) {
            strcpy(referer, p2);
          }
        }
        p = strtok(NULL, "\r\n");
      }

      if (uri[0] != '/') {
        strcpy(rescode, "400");
        strcpy(shortmsg, "Bad Request");
        strcpy(longmsg, "Bad URI");
        clienterror(connfd, ip, rescode, header, referer, u_agent, shortmsg, longmsg);
        close(connfd);
        continue;
      }
      querylen = strlen(uri+1);
      if (querylen > 500) {
        strcpy(result, "ERROR : Query sequence is too long. It must be 500 or less.");
      }
      else {
        strcpy(query, uri+1);
        for (i=0; i<querylen; i++) {
          query[i] = toupper(query[i]);
        }
        r = centroidfold(query, result, GAMMA);
        if (r) {
          strcpy(result, "ERROR : CentroidFold failed to execute this request.");
        }
      }
      strcat(body, result);
      strcat(body, "\n");
    }
    else if (!strcasecmp(method, "POST")) {
      reqbodylen = 0;
      strcpy(boundary, "\0");

      strcpy(tokwork,request);
      p = strtok(tokwork, "\r\n");
      if (p != NULL) {
        if (strlen(p) <= BUFSIZE) {
          strcpy(header, p);
        }
      }
      while(p != NULL) {
        if (!strncmp(p,"Content-Length", 14)) {
          if ((p2 = strpbrk(p, "123456789")) == NULL) {
            reqbodylen = 0;
          } else {
            reqbodylen = atoi(p2);
          }
        } 
        else if (!strncmp(p,"Content-Type", 12)) {
          p2 = strstr(p, "boundary=");
          if (p2 != NULL) {
            len = strlen(p2+9);
            if (strlen(p2+9) <= BUFSIZE) {
              strcpy(boundary, p2+9);
            }
            p2 = boundary + strlen(boundary) - 1;
            while (1) {
              if ((*p2 == '\r') || (*p2 == '\n')) { 
                *p2 = '\0';
                p2 --;
              } else {
                break;
              }
            }
          }
        }
        else if (!strncmp(p,"User-Agent", 10)) {
          p2 = p + 11;
          while (1) {
            if ((*p2 == ' ')) { 
              p2 ++;
            } else {
              break;
            }
          }
          len = strlen(p2);
          if (strlen(p2) <= BUFSIZE) {
            strcpy(u_agent, p2);
          }
        }
        else if (!strncmp(p,"Referer", 7)) {
          p2 = p + 8;
          while (1) {
            if ((*p2 == ' ')) { 
              p2 ++;
            } else {
              break;
            }
          }
          len = strlen(p2);
          if (strlen(p2) <= BUFSIZE) {
            strcpy(referer, p2);
          }
        }
        p = strtok(NULL, "\r\n");
      }
        
      if (reqbodylen == 0) {
        strcpy(rescode, "411");
        strcpy(shortmsg, "Length Required");
        strcpy(longmsg, "");
        clienterror(connfd, ip, rescode, header, referer, u_agent, shortmsg, longmsg);
        close(connfd);
        continue;
      }
      if (reqbodylen > CONTENTSIZE) {
        strcpy(body, "ERROR : Total size of query is too big. It must be 5K or less.");
        len = 0;
        while (1) {
          FD_ZERO(&mask);
          FD_ZERO(&readOK);
          FD_SET(connfd,&readOK);
          tv.tv_sec = 1;
          tv.tv_usec = 0;
          r = select(connfd+1, &readOK, NULL, NULL, &tv);
          if (FD_ISSET(connfd, &readOK)) {
            r = read(connfd,buf,BUFSIZE);
            len += r;
            buf[r] = '\0';
            if (len > MAXREADSIZE) {
              break;
            }
            if (r <= 0) {
              break;
            }
          } else {
            break;
          }
        }
        if (r < 0) {
          close(connfd);
          continue;
        }
      } else {
        p = strstr(request, "\n\r\n");
        if (p == NULL) {
          p = strstr(request, "\r\r\n");
        }
        p += 3;
        len = strlen(p);
        while (len < reqbodylen) {
          FD_ZERO(&mask);
          FD_ZERO(&readOK);
          FD_SET(connfd,&readOK);
          tv.tv_sec = TIMEOUT;
          tv.tv_usec = 0;
          r = select(connfd+1, &readOK, NULL, NULL, &tv);
          if (FD_ISSET(connfd, &readOK)) {
            r = read(connfd,buf,BUFSIZE);
            if (r <= 0) {
              break;
            }
            buf[r] = '\0';
            if (strlen(request)+r-1 > REQUESTSIZE) {
              len = -1;
              break;
            }
            strcat(request,buf);
          } else {
            break;
          }
          len = strlen(p);
        }
        if (r < 0) {
          close(connfd);
          continue;
        }
        if (len != reqbodylen) {
          strcpy(rescode, "411");
          strcpy(shortmsg, "Length Required");
          strcpy(longmsg, "");
          clienterror(connfd, ip, rescode, header, referer, u_agent, shortmsg, longmsg);
          close(connfd);
          continue;
        }

        if (strlen(boundary) != 0) {
          mp1 = strstr(p, boundary);
          if (mp1 != NULL) {
            mp2 = strstr(mp1+strlen(boundary), boundary);
            if (mp2 != NULL) {
              *mp2 = '\0';
              while (1) {
                mp2 --;
                if (*mp2 != '-') {
                  break;
                }
                *mp2 = '\0';
              }
            }
          }
        }

        mp1 = strchr(p+strlen(boundary), '>');
        while (1) {
          if (mp1 == NULL) {
            break;
          }
          mp2 = strchr(mp1, '\n');
          mp3 = strchr(mp1, '\r');
          if ((mp2 == NULL) && (mp3 == NULL)) {
            break;
          }
          if ((mp2 != NULL) && (mp3 != NULL)) {
            if (mp2 > mp3) {
              mp2 = mp3;
            }
          }
          if (mp2 == NULL) {
            mp2 = mp3;
          }
          mp3 = strchr(mp2, '>');
          if (mp3 == NULL) {
            mp3 = strchr(mp2, '\0');
          }

          startp = mp2 - p + 1;
          endp = startp + mp3 - mp2 - 1;
          strcpy(query, "\0");
          strcpy(result, "\0");
          querylen = 0;
          for (i=startp; i<endp; i++) {
            p[i] = toupper(p[i]);
            switch (p[i]) {
              case ' ':
              case '\t':
              case '\n':
              case '\r':
                break;
              default :
                if (querylen >= 500) {
                  strcpy(result, "ERROR : Query sequence is too long. It must be 500 or less.");
                }
                else {
                  query[querylen] = p[i];
                  querylen ++;
                }
            }
            if (strlen(result) != 0) {
              break;
            }
          }
          query[querylen] = '\0';
          if (strlen(result) == 0) {
            r = centroidfold(query, result, GAMMA);
            if (r) {
              strcpy(result, "ERROR : CentroidFold failed to execute this request.");
            }
          }
          strcat(body, result);
          strcat(body, "\n");
          mp1 = mp3;
          if (mp1 == '\0') {
            mp1 = NULL;
          }
        }
      }
    }
    else {
      strcpy(rescode, "501");
      strcpy(shortmsg, "Not Implemented");
      strcpy(longmsg, "Server dose not implement this method");
      clienterror(connfd, ip, rescode, header, referer, u_agent, shortmsg, longmsg);
      close(connfd);
      continue;
    }

    /* print response */
    t_byte = strlen(body);
    sprintf(buf, "HTTP/1.1 200 OK\nServer: CentroidFold Web Server\nContent-length: %d\nContent-type: text/plain\n\r\n", t_byte);
    r = write(connfd, buf, strlen(buf));
    if (r > 0) {
      //strcat(body,"\r\n");
      write(connfd, body, t_byte);
    }
    if (strcmp(SYSLOG, "none")) {
      tp = time(NULL);
      l_tp = localtime(&tp);
      if (strncmp(body, "ERROR", 5)) {
        sprintf(msg, "\t%s\t%s\t200\t%d\t%s\t%s\t%s\t%d", ip, header, t_byte, referer, u_agent, "OK", l_tp->tm_year+1900);
      } else {
        sprintf(msg, "\t%s\t%s\t200\t%d\t%s\t%s\t%s\t%d", ip, header, t_byte, referer, u_agent, body, l_tp->tm_year+1900);
      }
      syslog(LOG_INFO, msg); 
    }

    /* clean up */
    close(connfd);
  }

  /* close syslog */
  if (strcmp(SYSLOG, "none")) {
    closelog();
  }
}
